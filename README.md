# mTRFToolbox_simulations

The code in this repository is used to simulate the utility and test the effectiveness of the mTRF-Toolbox. To use this code, put it in a folder in the main directory containing the [mTRF-Toolbox][1]. The programs start with adding the `mtrf` path in the toolbox, which contains all the necessary functions.

These simulations were coded while writing an upcoming review on linear modeling for EEG and the mTRF-Toolbox:

> Crosse et al (in prep), "Linear modeling of neurophysiological responses to naturalistic stimuli: methodological considerations for applied research"

Specifically, Figure 3 of that paper is the output of the code in `Plot_TRFestimate_SNRnTr.m` which plots the results of the "high frequency" simulation (2-15 Hz) that is run in `Test_TRFestimate_SNRnTr_highfreq.m`. The code takes a while to run, so I have included a data file containing the results of the simulation plotted in Figure 3: `TRFestimate_w_SNRnTr_highfreq.mat`. Because the program takes a while to run, I ran it on the BlueHive server at the University of Rochester, and I have included the batch file to run it (`Test_TRFestimate_highfreq.sh`).

However, there were a number of different simulations I ran in addition to this one, and I have included those in this repository as well:

* `Test_TRFestimate_SNRnTr_delta.m`: Similar to the "high frequency" simulation, but I used a frequency range corresponding to delta (0.1 - 4 Hz). I have also included the related data file and batch file. To plot the results, I would change the name of the data file loaded in `Plot_TRFestimate_SNRnTr.m`.
* `Test_speechdataorig.m`: Inspired by the result that d-prime prediction accuracy could be used to validate the TRF model and its fit with the true model, I did forward envelope modeling of real EEG data collected during continuous speech listening (see the [Natural Speech dataset][2]) and looked at the model shapes as a function of d-prime prediction accuracy for one EEG channel (Fz). My hypothesis was that, if the models have good prediction accuracies, then they should look similar to each other if they represent real neural responses. I found that, for subjects with d-prime above 2, the models are similar, but for d-primes below 2 the models become more variable.

The results from these simulations / analyses were included in an earlier version of Figure 3, which I have included in this repository for reference.

Also:

* `LOO_nulltest_comparison.m`: LOO stands for "leave-one-out". Here, I was interested in examining how null distributions of prediction accuracies vary with the method of creating the null distribution. This was done in four different ways: 1) permute the testing trial, and compute null accuracy by pairing stimuli with true predictions (that is, calculated with the true model); 2) randomly circularly shift the stimulus in the testing trial, and compute null accuracy by pairing the shifted stimulus with the true prediction of the same trial; 3) permute all trials, but recalculate a new model on the permuted training trials before computing null accuracy (select the same lambda used for the true model); and 4) randomly circularly shift all trials, but recalculated a new model on the shifted training trials before computing null accuracy (using same lambda as true model). Our finding: it doesn't matter! The averages and variances of the null distributions are  comparable.

[1]: https://github.com/mickcrosse/mTRF-Toolbox
[2]:
