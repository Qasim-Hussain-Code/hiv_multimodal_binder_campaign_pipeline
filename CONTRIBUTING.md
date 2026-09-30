# Contributing

This repository has a house style and a few rules about evidence. A pull request that breaks one of them will be sent back, however good the science in it is.

## Writing rules

Every line of prose, every comment and every commit message should read as though a working scientist wrote it for other scientists.

1. No emojis anywhere: README, comments, `echo` output, commit messages.
2. No em dashes. Use a comma, a full stop, a colon or parentheses. Reshape the sentence instead of putting another character where the dash was going to be.
3. Banned phrases: "it is worth noting", "it is important to note", "in today's rapidly evolving", "plays a crucial role", "serves as a testament", "paving the way for", "in conclusion", "overall", "delve", "leverage" as a verb, "robust" as a compliment, "seamless", "comprehensive" as filler, "underscores", "showcases", "highlights" as a verb of significance.
4. No "not only X but also Y" and no "it is not X, it is Y". Both answer a misconception nobody had.
5. No trailing participial clause that tells the reader what a fact means. Write "Temsavir ranked 41st of 812." Stop there.
6. No bolded lead-in bullets that restate their own term.
7. Sentence case in headings.
8. Give numbers. A paragraph with no count, version, run time or threshold in it is probably filler.
9. Vary sentence length. Call arbitrary thresholds arbitrary. Say what did not work.
10. Comments explain reasoning, not syntax. Where a choice had an alternative, name it and say why it lost.

## Scores

A docking score is the output of an empirical function fitted to a training set. It is neither a free energy nor an affinity. Write the units and the function next to every score. Do not describe any computational score as an IC50 estimate, and do not write any sentence that could be read as a claim about treating or preventing infection in a person.

There is no pooled cross-modality score in this repository and a pull request adding one will be declined. The four modalities are ranked within themselves.

## Evidence

Every number in the README must trace to a named file in `results/` or `logs/`. If it does not, it does not go in the README. A pull request that adds a result must also add the script that produced it.

## The pre-registration

`results/preregistration.md` is committed once, on its own, with the message `add_preregistration`. It is never amended. A correction goes in a new file that names the original and its commit hash.

## Commits

Two or three lower case words joined with underscores. No prefixes, no scopes, no trailing punctuation. Examples: `add_ranking`, `fix_glycan_counts`.

## Code

Bash stages start with `set -euo pipefail`, print usage on `--help`, skip with a message when already complete, and register scratch directories for cleanup. Every stochastic step records its seed. Positive controls go through the same functions as every other candidate, and the code comment says so.

## Code of conduct

None is copied from elsewhere. Be courteous, criticise the work rather than the person, and assume the other author has run the thing.
