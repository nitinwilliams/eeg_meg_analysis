1. to install toolbox, add:
   a. file pipeline_masterscript.m to top of path
   b. folders HMM_functions, MPLS_functions to top of path

   this is done by lines 6-7 of pipeline_masterscript.m, if above file
   and folders are in current working directory.

2. to apply pipeline, run code in pipeline_masterscript.m - all other
   functions are called from there. 

3. if parallel computing available, adapt code in lines 11-13 of
   pipeline_masterscript.m to start parallel pool. if parallel 
   computing not available, change all 'parfor' loops to 'for' loops.

4. line 20 indicates names of files which point to input datasets.
   for each of these filenames, the variable 'bcidata_3D' contains the
   input data (as specified in line 34). change the filenames/variables
   to correspond to your input data. note that the input dataset should
   be in format: channels x samples x trials of both conditions

   it is assumed that there are 2 conditions, that they have same
   number of trials and that in the input dataset, all trials from one
   condition are followed by all trials from other condition.

5. lines 36-38 are specific to our BCI dataset and can be commented out.
   note that to run line 42 and onwards, the variable bcidata_3D needs 
   to contain input data in the format: 
   channels x samples x trials of both conditions

6. lines 42 to 55 perform baseline correction. 'bpts' specifies number of samples
   in pre-trial baseline window and 'tpts' specifies number of samples from stimulus
   onset to end of trial. bpts + tpts must be equal to the number of samples per trial,
   i.e. the size of second dimension of variable bcidata_3D

 
7. lines 59 to 65 is code for z-scoring each trial

8. lines 69 to 88 is code for z-scoring across trials

9. lines 112 to 148 segments dataset into a set of 'windows', which will also be the 
   input to the sparse-MVAR stage. setting parameter values correctly at this stage is
   especially important (see below):

   'srate' is the sampling frequency in Hz, 'M' is the number of channels. 'tpseg' is 
   the number of trials to be concatenated to produce a 'segment' for input to 
   sparse-MVAR stage. from each of these trials, only a small number of samples are 
   taken to produce the segment. 'seglen' is the total number of samples in a segment, 
   expressed in milliseconds, i.e. 'tpseg' multiplied by width of 1 trial in 
   milliseconds.

   in general, around 120 samples are required to estimate a sparse-MVAR model with
   28 signals and order 1. in the BCI data, we wanted each segment to be 50 ms wide.
   since sampling frequency is 240 Hz, this corresponded to 12 samples per segment
   (for each trial). since we needed 120 samples, we concatenated 10 trials. hence, 
   'seglen' is 500, i.e. 50 (milliseconds) x 10 trials.

10. lines 152-170 estimate sparsity level parameter for 10% of randomly selected 
    segments, for model orders from 1 to 6. this percentage is specified in pctval.

11. lines 174-201 determine optimal model order according to BCI, by estimating mode
    of model order estimated from 10% of randomly selected segments (given by pctval_2).

12. lines 205-221 estimate the sparse-MVAR model for each segment, with the selected
    sparsity level and model order. RMS value of sparse-MVAR coeffs. is taken.

13. for each segment, RMS value of s-MVAR coefficients are entered in vectorised form 
    in variable CMAT, which is a 2D matrix with dimensions:
    no. of segments x no. of s-MVAR coeffs. per segment

14. from line 239 onwards, variable CMAT is subject to clustering via k-means algorithm
    where k is varied according to the range previously specified.

15. from lines 243-251, the output of the clustering stage is re-ordered into the
    CONDMAT cell array, which contains two cells, i.e. one for each condition.
    each cell contains a 2D matrix wherein each 'segment' which was input to the 
    sparse-MVAR stage is replaced by a 'symbol' which is the output of the clustering.

16. lines 259-281 perform Markov Modelling on 'symbols' from both conditions
    separately, after initialising the emission matrix to identity in both cases.

17. from lines 285-338, the performance measure, i.e. model distance is also
    assessed for statistical significance, to determine if Markov Models estimated from both
    conditions are indeed different. null distribution for 'model distance' is also 
    determined using a permutation-based procedure. model distance for each subject, 
    for each value of 'k' is in variable MD

19. lines 340-349 generate p-values of MD for each subject for each value of 'k'
    p-values for MD are in 'MD_pvals'

20. hope this helps! please e-mail nitinwilliams@gmail.com if you have any questions!
         
   