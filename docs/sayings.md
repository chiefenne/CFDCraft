# Sayings

Small words from the shop floor of CFD, collected for the community. Add yours with a pull request.

## Residuen-gaffen

*Austrian German, "Residuen" = residuals, "gaffen" = dialect for staring.* The pose of the person who sits in front of the screen and watches the residuals of a run go by, number after number. A heavy transient run is the best occasion for it, and the reason why the log shows the residuals of **every outer iteration** under the box of each logged time step (see [the log](gallery/runs.md#log)): it gives the one who stares something to see.

A secretary once asked the people at the screens whether they typed these numbers in on their own. Nobody did, the solver wrote them. Now the log is made for people who like to look at them, and the machine files (`run.jsonl`) for the ones who would rather not.

## Deemed prudent

*English, the voice of a report.* "It is deemed prudent to investigate the nature of heat transfer in more depth." Written by a native speaker of English who knew nothing about heat transfer, and every word of it is right: the grammar, the tone, the careful caution. It says nothing, and it says it perfectly. It also fits any report on any topic: put "turbulence" or "the pressure drop" in place of "heat transfer" and it still reads well.

Foehncraft takes the advice anyway: energy, buoyancy and conjugate heat transfer are stage 7.

## Residuals got stuck

*English, from a solver developer.* Said of a run whose residuals stopped falling and stayed flat, as if they were a car in the snow. The words come from someone who knows exactly what is behind it, and still they make the residuals sound like the guilty party. They are only the messenger: a flat line points to the mesh, the boundary conditions, the relaxation or the physics itself, and the residuals of every outer iteration in the log are there to find where it began.

## Better to know what's not happening

*English, from an R&D manager.* "It is better to know what's not happening, than not knowing what's happening." Said about a computation, by a manager who did not trust simulation yet and was sure he could check the CFD code with Excel. A few cells for the mass flow, a few for the pressure drop, and the spreadsheet would show whether the solver told the truth.

The sentence turns on itself twice and still lands. Read it slowly and it is a fair demand on any solver: a run should say plainly what it does not do, and never leave you guessing what it does. That is why Foehncraft writes the balances of every run into the ledger ([gallery 19](gallery/runs.md#health)), checks the boundary conditions before it starts and keeps a [verification manual](verification/manual.md): so that nobody needs a spreadsheet to find out what is happening, and what is not.
