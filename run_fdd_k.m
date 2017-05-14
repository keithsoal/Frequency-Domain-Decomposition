% Call and run fdd_k.m
clear;close all;clc;

load('measurement1.mat');

% *********************************************
%               Measurement parameters 
% *********************************************
fs =4096;                                   % Sample frequency
nfft=4096;                                  % Number of fft points
delta_freq = fs/nfft;                       % Frequency resolution
w = hann(nfft);                             % Create window

[Frq,phi] = fdd_k(y,fs,nfft,w);