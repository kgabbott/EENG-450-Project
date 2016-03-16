# EENG 450 Project
This is the code for my EENG 450 (Applied Digital Signal Processing Final Project). The purpose of this project is to replicate the 
speech encoding perfromed by telephone communication. In this process, speech is converted into coeffitients of a linerar model, which
is then used to sythesize the speech at the reciever.

The system in this project is able to take in a speech segment, encode it, and play back a sysnthesized version. The encoding proess 
involves taking a 50ms segment, and using the autocorrelation and spectrum to determine the pitch frequency and locate the poles/zeros 
for a 6th order ARMA model. This model is then driven by a pitch impulse train (if voiced speech) or random noise (if unvoiced speech). 
