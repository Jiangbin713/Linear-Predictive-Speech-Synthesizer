%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2018 assignment of Speech & Audio Processing & Recognition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BIN JIANG 6519680
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Important Variants:
%s1 : Female voice source data
%s2 : Male voice source data
%a1 : LPC coefficient of Female speech
%a2 : LPC coefficient of Male speech
%formant1 : Formant of Female speech
%formant2 : Formant of Male speech
%f_F_average1: Fundamental frequency of female speech
%f_F_average2: Fundamental frequency of male speech
%synthesize_vowel_1 : Synthesized vowel of female speech
%synthesize_vowel_2 : Synthesized vowel of male speech
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Import commands:
%sound(s1,fs1) play female voice 'head'
%sound(s2,fs2) play male voice 'head'
%sound(synthesize_vowel_1) play synthesized vowel of female voice
%sound(synthesize_vowel_2) play synthesized vowel of male voice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;

% %---------------------------------------------------------------------
% %---------------------------------------------------------------------
% % Female Voice 'head'
% %-------------------------------------------------------------------------
% %-------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Load female voice :'head'
[s1,fs1] = audioread('head_f.wav'); %sample data & sample rate
%sound(s1); %play the sound

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Plot the original speech waveform
figure(1)
subplot(2,1,1)
tlabel_1 = (0:1:length(s1)-1).*(1/fs1);
plot(tlabel_1,s1);  %plot original Female sound
title('Female voice waveform');xlabel('Time (s)'); ylabel('Magtitude');%label
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Plot the speech amplitude spectrum
n1 = length(s1);
sfft1 = fft(s1); %samplefft   fft length = sample length
sfft1 = sfft1((1:floor(n1/2))); % because of the symmtery
fftf1 = (1:n1/2)*(fs1/n1); % because of the symmtery
mag1 = abs(sfft1);
subplot(2,1,2)
plot(fftf1,20*log10(mag1)) %frequency domain graph
title('Spectrum of the female speech Segment'); xlabel('Frequency (Hz)') ; ylabel('Magtitude(dB)'); %label
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Calculate the LPC coefficient and Plot the filter resonse
a1 = lpc(s1,128); %LPC coefficient 
a11 = lpc(s1,20);
a12 = lpc(s1,6);
[H1,F1] = freqz(1, a1, floor(n1/2), fs1); 
[H11,F11]=freqz(1, a11, 1:floor(n1/2), fs1);
[H12,F12]=freqz(1, a12, floor(n1/2), fs1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Plot the filter response
plot(F1,20*log10(abs(H1))) 
plot(F11,20*log10(abs(H11)), 'k')
plot(F12,20*log10(abs(H12)), 'y')
legend('Amplituide spectrum', 'LPCspecturm, p=128','LPCspecturm, p=20','LPCspecturm, p=6');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Spectrogram of Female Voice
figure(2)
subplot(2,1,1)
spectrogram(s1,400,[],[],'yaxis')
title('Spectrogram of Female Voice');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Formant Frequency Estimation
%s1window = s1.*hamming(n1);
%s1window = filter(1, [1 0.63] ,s1window);
a11 = lpc(s1,20);
rts1 = roots(a11);
rts1 = rts1(imag(rts1) >=0);
angz1 = atan2(imag(rts1), real(rts1));
[frqs1, indices1] = sort(angz1.* (fs1/(2*pi)));


num = 1;
for c = 1 : length(frqs1)
    if (frqs1(c) > 0)
        formant1(num) = frqs1(c);
        num = num+1;
    end
end
% %Print the first three formant frequency
fprintf('\n %s :\n','The first three formant frequency of female voice')
fprintf('%3f Hz \n',formant1(1:3))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  Fundamental Frequency Estimation
fundamental_F1 = pitch(s1,fs1);
[F_amount1,y] = size(fundamental_F1); 
f_F_average1 = sum(fundamental_F1)/F_amount1; %mean Fundamental Frequency
fprintf('\n %s :\n','The mean fundamental frequency of female voice')
fprintf('%3f Hz \n',f_F_average1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Synthesized vowels
p1 = [];
n = 1;
for t = 1: 1 :fs1   % generate impluse train around 
    if t == n*106  %106 is calculated by hand
        p1(t) = 0.3;
        n = n+1;
    else
        p1(t) = 0;
    end
end
t = 1 : 1: n1;  
tp_1 = (0:1:fs1-1).*(1/fs1);
figure(3)
subplot(2,1,1)
plot(tp_1,p1)  %plot the train impluse of male voice
axis([0 inf 0 1])
title('Female train impluse'); xlabel ('Time (s)'); ylabel('Magtigude')
hold on
synthesize_vowel_1 = filter(a1,1,p1); %generate the synthesized male vowel
%sound(synthesize_vowel_1);

figure(1) % plot the synthesized vowel together with source data
subplot(2,1,1)
plot(tp_1(1:n1),synthesize_vowel_1(1:n1),'--')
legend('Original female voice', 'Synthesize female vowel')

%audiowrite('synthesized_female_head.wav' ,synthesize_vowel_1,24000);





% %---------------------------------------------------------------------
% %---------------------------------------------------------------------
% %Male Voice 'head'
% %---------------------------------------------------------------------
% %---------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Load male voice :'head'
[s2,fs2] = audioread('head_m.wav'); %sample data & sample rate
%sound(s2,fs2); %plot original Male sound

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Plot the original speech waveform
n2 = length(s2');
tlabel_2 = (0:1:length(s2)-1).*(1/fs2);
figure(4)
subplot(2,1,1)
plot(tlabel_2,s2); %plot original Male sound 
title('Male voice waveform ');xlabel('Time (s)'); ylabel('Magtitude'); %label
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Plot the speech amplitude spectrum
sfft2 = fft(s2,n2); %samplefft   fft length = sample length
sfft2 = sfft2(floor(1:(n2/2))); % because of the symmtery
fftf2 = (1:floor(n2/2))*(fs2/n2); % because of the symmtery
mag2 = abs(sfft2);
subplot(2,1,2)
plot (fftf2,20*log10(mag2)) %frequency domain graph
title('Spectrum of the male speech Segment'); xlabel('Frequency (Hz)') ; ylabel('Magtitude(dB)');%label
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Calculate the LPC coefficient and Plot the filter resonse
a2 = lpc(s2,128); %LPC coefficient 
a21 = lpc(s2,20);
a22 = lpc(s2,6);
[H2,F2] = freqz(1, a2, floor(1:(n2/2)), fs2); 
[H21,F21]=freqz(1, a21, floor(1:(n2/2)), fs2);
[H22,F22]=freqz(1, a22, floor(1:(n2/2)), fs2);
 
plot(F2,20*log10(abs(H2))) %Plot the filter response
plot(F21,20*log10(abs(H21)), 'k')
plot(F22,20*log10(abs(H22)), 'y')
legend('Amplituide spectrum', 'LPCspecturm, p=128','LPCspecturm, p=20','LPCspecturm, p=6');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Spectrogram Male Voice
figure(2)
subplot(2,1,2)
spectrogram(s2,400,[],[],'yaxis')
title('Spectrogram of Male Voice');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Formant Frequency Estimation
%s2window = s2.*hamming(n2);
%s2window = filter(1, [1 0.63] ,s2window);
a22 = lpc(s2, 20);
rts2 = roots(a22);
rts2 = rts2(imag(rts2) >=0);
angz2 = atan2(imag(rts2), real(rts2));
[frqs2, indices2] = sort(angz2.* (fs2/(2*pi)));


num = 1;
for c = 1 : length(frqs2)
    if (frqs2(c) > 0)
        formant2(num) = frqs2(c);
        num = num+1;
    end
end
formant2(1:3);% the first three formant frequency
fprintf('\n %s :\n','The first three formant frequency of male voice')
fprintf('%3f Hz \n',formant2(1:3))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  Fundamental Frequency Estimation 
fundamental_F2 = pitch(s2,fs2);
[F_amount2,y] = size(fundamental_F2); 
f_F_average2 = sum(fundamental_F2)/F_amount2;
fprintf('\n %s :\n','The mean fundamental frequency of male voice')
fprintf('%3f Hz \n',f_F_average2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Synthesized vowels

p2 = [];
n = 1;
for t = 1: 1 :fs2   % generate impluse train around 1s
    if t == n*186  %186 is calculated by hand
        p2(t) = 1;
        n = n+1;
    else
        p2(t) = 0;
    end
end
t = 1 : 1: n2;  
tp_2 = (0:1:fs2-1).*(1/fs2);
figure(3)
subplot(2,1,2)
plot(tp_2,p2) 
title('Male train impluse'); xlabel ('Time (s)'); ylabel('Magtigude')
synthesize_vowel_2 = filter(1,a2,p2);  %generate the synthesized female vowel
%sound(synthesize_vowel_2);
figure(4)
subplot(2,1,1)
plot(tlabel_2(1:n2),synthesize_vowel_2(1:n2),'--') % plot the synthesized vowel together with source data
legend('Original male voice', 'Synthesize male vowel')

%audiowrite('synthesized_male_head.wav' ,synthesize_vowel_2,24000);

