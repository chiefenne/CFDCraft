# Sayings

Small words from the shop floor of CFD, collected for the community. Add yours with a pull request.

## Residuen-gaffen

*Austrian German, "Residuen" = residuals, "gaffen" = dialect for staring.* The pose of the person who sits in front of the screen and watches the residuals of a run go by, number after number. A heavy transient run is the best occasion for it, and the reason why the log shows the residuals of **every outer iteration** under the box of each logged time step (see [the log](gallery/runs.md#log)): it gives the one who stares something to see.

A secretary once asked the people at the screens whether they typed these numbers in on their own. Nobody did, the solver wrote them. Now the log is made for people who like to look at them, and the machine files (`run.jsonl`) for the ones who would rather not.
