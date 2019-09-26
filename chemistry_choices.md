The aim of this file is to document choices that handled bug correction and why, concerning precise chemoinformatics choices.

# Chemical rule application.

When a rule applies to a substrate and after standardisation produces this substrate again (S -> S + I), the rule is deleted as this is not biological.
This is corrected at the compound stage.

# Compound equality: either main layer or full inchikey
Choices: usually less stringent for the chassis.

# Moves generating duplicate compounds:
- Only unique compounds are conserved.
- Logs will say it is merged (and conserve the number of compounds in stoechiometry dictionnary)


# Moves generating the same compounds:

Keep the one with the higher score. In practice, also keeping the synonyms (transformation IDs) of the other moves generating the same compounds.

# History of the state.

Keeping a history of visited compounds (excluding the organism's compounds).
- Refuse moves that generate compounds present in the history to avoid loops.

# Refusing rules that produce a different number of compounds as the original.

This can happen when the rule that learned on 2 molecules, which subgroups appear in a much bigger molecule when doing the retrosynthetic search.
It is unrealistic to expect an enzyme to work this way.
