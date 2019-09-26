r"""
Apply one/many rules on one/many substrates.

Thomas Duigou, INRA, 2018
"""

import os
import sys
import gzip
import json
import time
import signal
import logging
import argparse
import multiprocessing as mp

from rdkit import Chem
from rdkit.Chem import AllChem

from utilities.reactor.Core import RuleBurnerCore
from utilities.reactor.Utils import standardize_chemical, standardize_results, handle_results, ChemConversionError


class RuleConversionError(Exception):
    """Raised when something went wrong during SMARTS conversion to RDKit rxn object."""

    def __init__(self, msg):
        self._msg = msg

    def __str__(self):
        return "RULE-CONVERSION-ERROR: {}".format(self._msg)


def worker_match(kwargs):
    """Check if a chemical can be fired by a rule according to left side."""
    w = RuleBurnerCore(**kwargs)
    return w.match()


def worker_fire(kwargs):
    """Apply a reaction a rule on a chemical."""
    r = RuleBurnerCore(**kwargs)
    return r.fire()


def kill(pool):
    """Send SIGTERMs to kill all processes belonging to pool.

    Will not work on Windows OS.
    """
    # stop repopulating new child
    pool._state = mp.pool.TERMINATE
    pool._worker_handler._state = mp.pool.TERMINATE
    for p in pool._pool:
        os.kill(p.pid, signal.SIGKILL)
    # .is_alive() will reap dead process
    while any(p.is_alive() for p in pool._pool):
        pass
    # Get-lucky workaround: force releasing lock
    try:
        pool._inqueue._rlock.release()
    except ValueError as e:
        logging.error(e)
    pool.terminate()


class RuleBurner(object):
    """Apply any number of rules on any number of compounds."""

    def __init__(
            self, rsmarts_list, inchi_list, rid_list=None,  cid_list=None,
            match_timeout=1, fire_timeout=1, ofile=None, compress=False):
        """Setting up everything needed for behavor decisions and firing rules.

        :param  rsmarts_list:   list of reaction rule SMARTS
        :param  inchi_list:     list of inchis
        :param  rid_list:       list of reaction rule IDs
        :param  cid_list:       list of chemical IDs
        :param  match_timeout:  int, timeout execution for compound pre-matching
        :param  fire_timeout:   int, timeout execution for rule firing
        :param  ofile:          str, Output file to store results
        """

        # Internal settings
        self._INDENT_JSON = True
        self._TRY_MATCH = False

        # Input
        self._rsmarts_list = rsmarts_list
        self._rid_list = rid_list
        self._inchi_list = inchi_list
        self._cid_list = cid_list

        # Settings
        self._try_match = self._TRY_MATCH  # TODO: add option for that
        self._match_timeout = match_timeout
        self._fire_timeout = fire_timeout

        # Check for consistency between depictions and IDs
        try:
            if self._rid_list:
                assert len(self._rsmarts_list) == len(self._rid_list)
        except AssertionError as e:
            logging.warning("ID and depiction rule lists have different size, compound IDs will be ignored")

        try:
            if self._cid_list:
                assert len(self._inchi_list) == len(self._cid_list)
        except AssertionError as e:
            logging.warning("ID and depiction compounds lists have different size, rule IDs will be ignored")

        # Output
        self._json = list()
        self._compress = compress
        if not ofile:
            self._ofile = None
        else:
            pdir = os.path.abspath(os.path.dirname(ofile))
            os.makedirs(pdir, exist_ok=True)
            self._ofile = os.path.abspath(ofile)

        # A place to swim
        self._pool = None

    def _run_with_timeout(self, worker, kwargs, timeout=5):
        """Generic wrapper making use of multiprocessing to garantee effective timeout.

        :param  worker:     function to be called
        :param  kwargs:     dictionnary of args to be passed to the called function
        :param  timeout:    int, timeout
        :returns            depends of the worker function
        """
        if self._pool is None:
            self._pool = mp.Pool(processes=1)
        try:
            start_time = time.time()
            res = self._pool.apply_async(worker, kwds={'kwargs': kwargs})
            ans = res.get(timeout=timeout)
            end_time = time.time()
            exec_time = round(end_time - start_time, 4)
        except mp.TimeoutError as e:
            kill(self._pool)
            self._pool = mp.Pool(processes=1)
            raise e
        except Exception as e:
            raise e
        return ans, exec_time

    def _jsonify(self, rsmarts, inchi, rid=None, cid=None,
                 has_match=None, match_timed_out=None,
                 match_exec_time=None, match_error=None,
                 fire_timed_out=None,
                 fire_exec_time=None, fire_error=None,
                 inchikeys_list=None, inchis_list=None, smiles_list=None):
        """Return the results as a JSON string.

        :param      rsmarts:            str, reaction rule string depiction
        :param      inchi:              str, substrate string depiction
        :param      rid:                str, reaction rule ID
        :param      cid:                str, substrate ID
        :param      has_match:          bool or None, True if there is match
        :param      match_timed_out:    bool, True if timeout reached
        :param      match_exec_time:    int, execution time for matching
        :param      match_error:        str, error message if any, else None
        :param      fire_timed_out:     bool, True if timeout reached
        :param      fire_exec_time:     bool, execution time for firing
        :param      fire_error:         str, error message if any, else None
        :param      inchikeys_list:     list of list, Inchikeys of products
        :param      inchis_list:        list of list, Inchis of products
        :param      smiles_list:        list of list, SMILES of products
        :returns    json_string:        JSON string
        """
        # General info
        data = {
            'rule_id': rid,
            # 'rule_smarts': rsmarts,
            'substrate_id': cid,
            # 'substrate_inchi': inchi,
        }
        # Match info
        if self._try_match:
            data['match'] = has_match
            data['match_timed_out'] = match_timed_out
            data['match_exec_time'] = match_exec_time
            if match_error is not None:
                data['match_error'] = match_error
        # Fire info
        data['fire_timed_out'] = fire_timed_out
        data['fire_exec_time'] = fire_exec_time
        if fire_error is not None:
            data['fire_error'] = fire_error
        if (inchikeys_list is not None) and (len(inchikeys_list) > 0):
            data['product_inchikeys'] = inchikeys_list
            data['product_inchis'] = inchis_list
            data['product_smiles'] = smiles_list

        return json.dumps(obj=data, indent=self._INDENT_JSON)

    def write_json(self):
        """Write the JSON string."""
        # Handling file handler
        if self._ofile:
            if self._compress:
                ofh = gzip.open(self._ofile, 'wb', compresslevel=9)
            else:
                ofh = open(self._ofile, 'w')
        else:
            ofh = sys.stdout
        # Big string
        content = '[\n' + ','.join(self._json) + '\n]' + '\n'
        # Handling compression
        if self._ofile and self._compress:
            ofh.write(content.encode())
        else:
            ofh.write(content)
        ofh.close()

    def compute(self):
        """Rules under fire."""
        for rindex, rsmarts in enumerate(self._rsmarts_list):
            # Extract corresponding reaction rule ID if any
            if self._rid_list:
                rid = self._rid_list[rindex]
            else:
                rid = None
            # Get RDKit reaction object
            try:
                rd_rule = AllChem.ReactionFromSmarts(rsmarts)
                rd_rule.Initialize()
            except Exception as e:
                raise RuleConversionError(e) from e

            for cindex, inchi in enumerate(self._inchi_list):
                # Extract corresponding substrate ID if any
                if self._cid_list:
                    cid = self._cid_list[cindex]
                else:
                    cid = None
                # Get standardized RDKit mol
                try:
                    # rd_mol = Chem.MolFromSmiles(csmiles, sanitize=False)  # Important: Sanitize = False
                    rd_mol = Chem.MolFromInchi(inchi, sanitize=False)  # Important: Sanitize = False
                    rd_mol = standardize_chemical(rd_mol, add_hs=False, rm_stereo=True, heavy=True)
                except Exception as e:
                    raise ChemConversionError(e) from e
                # General args to used for both matching and firing
                kwargs = {
                        'rd_rule': rd_rule,
                        'rd_mol': rd_mol
                        }
                # Matching
                has_match = None
                match_exec_time = None
                match_timed_out = None
                match_error = None
                if self._try_match:
                    try:
                        has_match, match_exec_time = self._run_with_timeout(
                                worker=worker_match, kwargs=kwargs,
                                timeout=self._match_timeout
                                )
                        match_timed_out = False
                        match_error = None
                    except mp.TimeoutError as e:
                        has_match = None
                        match_exec_time = None
                        match_timed_out = True
                        match_error = str(e)
                    except Exception as e:
                        has_match = None
                        match_exec_time = None
                        match_timed_out = False
                        match_error = str(e)
                # Firing
                try:
                    ans, fire_exec_time = self._run_with_timeout(
                            worker=worker_fire, kwargs=kwargs,
                            timeout=self._fire_timeout
                            )
                    rdmols, failed = standardize_results(ans, add_hs=False, rm_stereo=True)
                    inchikeys, inchis, smiles = handle_results(rdmols)
                    fire_timed_out = False
                    fire_error = None
                except ChemConversionError as e:
                    inchikeys = None
                    inchis = None
                    smiles = None
                    fire_timed_out = False
                    fire_error = str(e)
                    logging.warning(e)
                except mp.TimeoutError as e:
                    fire_exec_time = None
                    inchikeys = None
                    inchis = None
                    smiles = None
                    fire_timed_out = True
                    fire_error = str(e)
                    logging.error('TIMEOUT: cid={}, rid={}'.format(cid, rid))
                    logging.error('TIMEOUT: original error={}'.format(e))
                except Exception as e:
                    fire_exec_time = None
                    inchikeys = None
                    inchis = None
                    smiles = None
                    fire_timed_out = False
                    fire_error = str(e)
                    logging.warning(e)
                # JSONify and store
                json_str = self._jsonify(
                        rsmarts=rsmarts,
                        inchi=inchi,
                        rid=rid,
                        cid=cid,
                        has_match=has_match,
                        match_timed_out=match_timed_out,
                        match_exec_time=match_exec_time,
                        match_error=match_error,
                        fire_timed_out=fire_timed_out,
                        fire_exec_time=fire_exec_time,
                        fire_error=fire_error,
                        inchikeys_list=inchikeys,
                        inchis_list=inchis,
                        smiles_list=smiles
                        )
                self._json.append(json_str)


def __cli():
    """Command line interface."""

    help = "Apply rules on chemicals."

    def inline_mode(args):
        """Execution mode to be used when a single rule and a single chemical
        are provided through CLI.
        """
        r = RuleBurner(
                rsmarts_list=[args.rsmarts], inchi_list=[args.inchi],
                rid_list=[args.rid], cid_list=[args.cid],
                match_timeout=args.match_timeout, fire_timeout=args.fire_timeout,
                ofile=args.ofile, compress=args.compress
                )
        r.compute()
        r.write_json()

    def infile_mode(args):
        """Execution mode to be used when rules and chemicals are provided
        in CSV files.
        """

        rsmarts_list = list()
        rids_list = list()

        inchi_list = list()
        cids_list = list()

        import csv

        with open(args.rfile, 'r') as ifh:
            reader = csv.DictReader(ifh, delimiter='\t')
            for row in reader:
                rsmarts_list.append(row['rule_smarts'].strip())
                rids_list.append(row['rule_id'].strip())

        with open(args.cfile, 'r') as ifh:
            reader = csv.DictReader(ifh, delimiter='\t')
            for row in reader:
                inchi_list.append(row['inchi'].strip())
                cids_list.append(row['chem_id'].strip())

        r = RuleBurner(
                rsmarts_list=rsmarts_list, inchi_list=inchi_list,
                rid_list=rids_list, cid_list=cids_list,
                match_timeout=args.match_timeout, fire_timeout=args.fire_timeout,
                ofile=args.ofile, compress=args.compress
                )
        r.compute()
        r.write_json()

    parser = argparse.ArgumentParser(description=help)
    parser.add_argument('--match_timeout', help='Rule matching timeout. Default: 1.', default=1, type=int)
    parser.add_argument('--fire_timeout', help='Rule furing timeout. Default: 1.', default=1, type=int)
    parser.add_argument('--ofile', help='Output file to store results. Default to STDOUT if none provided')
    parser.add_argument('--compress', action='store_true', help='Enable gzip compression (only when output to file).')
    subparsers = parser.add_subparsers(help='Input mode')

    parser_inline = subparsers.add_parser('inline', help='Get inputs from command line')
    parser_inline.set_defaults(func=inline_mode)
    parser_inline.add_argument('--rsmarts', help='Reaction rule SMARTS', required=True)
    # parser_inline.add_argument('--csmiles', help='Chemical SMILES depiction', required=True)
    parser_inline.add_argument('--inchi', help='Chemical inchi depiction', required=True)
    parser_inline.add_argument('--rid', help='Reaction rule ID, optional')
    parser_inline.add_argument('--cid', help='Chemical ID, optional')

    parser_file = subparsers.add_parser('infile', help='Get inputs from files')
    parser_file.set_defaults(func=infile_mode)
    parser_file.add_argument(
            '--rfile', required=True, help=' '.join([
                    'Reaction rule file.',
                    'Tab separated columns.',
                    'One reaction rule per line.',
                    'Mandatory column: rule_smarts.',
                    'Optional column: rule_id.',
                    'Other columns will be ignored.'
                    ])
            )
    parser_file.add_argument(
            '--cfile', required=True, help=' '.join([
                    'Chemical file.',
                    'Tab separated columns.',
                    'One chemical per line.',
                    'Mandatory column: inchi.',
                    'Optional column: chem_id.',
                    'Other columns will be ignored.'
                    ])
            )

    # Logging
    logging.basicConfig(
            stream=sys.stderr, level=logging.INFO,
            datefmt='%d/%m/%Y %H:%M:%S',
            format='%(asctime)s -- %(levelname)s -- %(message)s'
            )

    # Execute right mode
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    __cli()
