% Kevin Abbott
% EE450: Final Project

function Project
global seg_size hamm order pitch_modulation

order = 6;
seg_size = 400;

pitch_cutoff = 20;
energy_cutoff = .001;

isrt = 1;

pitch_modulation = 2;

i=0:seg_size-1;
hamm = 0.54 - 0.46*cos(2*pi*i/(seg_size-1));

figure(1), clf
figure(2), clf

exit = 0;
while exit == 0
    choice = menu('Project','Speech Acquisition','Speech Analysis','Speech synthesis','Exit');
    switch choice
        case 1
%           Record and play sample
            signal = record();
            toFile(signal, 'SpeechOrig.txt')
            soundsc(signal);
            
            figure(1), clf, title('Original Sample')
            plot(signal)
        case 2 
%           Encode speech sample
            x = fromFile('SpeechOrig.txt');
            soundsc(x);

            isrt = find(x>.2, 1);
            
            figure(1), subplot(211)
            plot(x)
            axis off, title('Original Speech')
            
            figure(2), subplot(211)
            plot(x(isrt:min(16000, isrt+1200)))
            axis off, title('Segment of Original Speech')
            
            fid = fopen('SpeechCode.txt','w');
%           Encode each 50ms segment of speech
            for i = 1:length(x)/seg_size
                segment = x((seg_size*(i-1))+1:seg_size*i);
                energy = (segment'*segment);
                if energy <= energy_cutoff
                   fprintf(fid, '0\n'); 
                else
                    energy = round(energy*20);
                    [pitch, rs, thetas] = process(segment);  
                    if pitch < pitch_cutoff
                        pitch = 0;
                    end
                    thetas = floor(thetas*255/pi);
                    rs = floor(rs*255);
                    fprintf(fid, '%i %i', energy, pitch);
                    for j = 1:3
                        fprintf(fid, ' %i %i', rs(j), thetas(j));
                    end
                    fprintf(fid, '\n');
                                        
                end
            end
            fclose(fid);
            
        case 3
%          Reconstruct encoded speech 
           pitch = 0;
           rs = [];
           thetas = [];
           
           synth = [];
           synth_p_increased = [];
           synth_p_decreased = [];
           
%          Get encoded values for each segment
           fid = fopen('SpeechCode.txt','r');
           line = fgetl(fid);
           while ischar(line)
               params = str2num(char(strsplit(line)));
               energy = params(1)/20;
               if energy ~= 0
                   pitch = params(2);
                   rs = [params(3) params(5) params(7)];
                   thetas = [params(4) params(6) params(8)];
                   rs = rs/255;
                   thetas = (thetas/255)*pi;
               end
%              Synthesize the speech at original and modified pitches
               [seg, seg_p_increased, seg_p_decresed] = synthesize(energy, pitch, rs, thetas);
               synth = [synth seg];
               synth_p_increased = [synth_p_increased seg_p_increased];
               synth_p_decreased = [synth_p_decreased seg_p_decresed];
               line = fgetl(fid);
           end
           fclose(fid);
           
%          Plot and play the synthesized speech at original, lower, and
%          higher pitches
           figure(1), subplot(212)
           plot(synth)
           axis off, title('Synthesized Speech')
           figure(2), subplot(212)
           plot(synth(isrt:min(16000, isrt+1200)))
           axis off, title('Segment of Synthesized Speech')
           sound(synth);
           pause(2)
           
           figure(1), subplot(212)
           plot(synth_p_increased) 
           axis off, title(sprintf('Synthesized Speech: Pitch Period*%d', pitch_modulation))
           figure(2), subplot(212)
           plot(synth_p_increased(isrt:min(16000, isrt+1200)))
           axis off, title(sprintf('Segment of Synthesized Speech: Pitch Period*%d', pitch_modulation))
           sound(synth_p_increased);
           pause(2)
           
           figure(1), subplot(212)
           plot(synth_p_decreased) 
           axis off, title(sprintf('Synthesized Speech: Pitch Period/%d', pitch_modulation))
           figure(2), subplot(212)
           plot(synth_p_decreased(isrt:min(16000, isrt+1200)))
           axis off, title(sprintf('Segment of Synthesized Speech: Pitch Period/%d', pitch_modulation))
           sound(synth_p_decreased);
        otherwise
            exit = 1;
            
    end
end  % while
end % function main


function [pitch, rs, thetas] = process(segment)
global hamm order

x_ham = hamm'.*segment;

% Calculate cepstrum from sepctral estimation
S = abs(fft(x_ham)).^2;
cepstrum = ifft(10*log10(S));
cepstrum = cepstrum(1:ceil(length(cepstrum)/2));
[~, pitch] = max(cepstrum(10:length(cepstrum)));
pitch = pitch + 9;

Rx = Acor(x_ham, order);
%Create the toeplitz matrix
R = toeplitz(Rx(1:order));
G = Rx(2:order+1)';
%Calculate the impulse response
H = R\G;
%Find the roots of the system
C = [1 -H'];

%   Find and the roots of the linear prediction filter
rts = roots(C);
thetas = [];
rs = [];
zero_theta = [];

for j = 1:order
    rt = rts(j,:);
    theta = atan2(imag(rt),real(rt));
    r = abs(rt);
    % If roots unstable, bring inside unit circle
    if r > 1
       r = 1/r;
    end
    
    if theta > 0
        thetas = [thetas theta];  
        rs = [rs r];
    elseif theta == 0
        zero_theta = [zero_theta r];
    end
end
if not(isempty(zero_theta))
    rs = [rs max(zero_theta)];
    thetas = [thetas 0];
end
num_poles = min(3, length(rs));
rs = rs(1:num_poles);
thetas = thetas(1:num_poles);
end

function [seg, seg_p_increased, seg_p_decresed] = synthesize(energy, pitch, rs, thetas)
global seg_size pitch_modulation

if energy == 0
%     If energy zero, just produce segment of silence
    seg = zeros(1, seg_size);
    seg_p_increased = seg;
    seg_p_decresed = seg;
else
    % Calculate the AR coeffitients
    AR_coeff = zeros(3,2);
    for i = 1:3
        AR_coeff(i,:) = [(2*rs(i)*cos(thetas(i))) -(rs(i))^2];
    end
    
    if pitch == 0
        % Synthesize the segment of unvoiced speech
        unvoiced = rand(1,seg_size);
        for i = 1:length(AR_coeff)
            unvoiced = AR(unvoiced, AR_coeff(i, :));
        end
        seg = unvoiced;
        seg_p_increased = unvoiced;
        seg_p_decresed = unvoiced;
    else
       % Synthesize the segment of voiced speech
       seg = zeros(1, seg_size); 
       for i = 0:round(seg_size/pitch)-1
           seg(i*pitch+1) = 1;
       end
       for i = 1:length(AR_coeff)
           seg = AR(seg, AR_coeff(i, :));
       end
       
       % Synthesize the segment of voiced speech with higher pitch period
       increased_p = round(pitch*pitch_modulation);
       seg_p_increased = zeros(1, seg_size); 
       for i = 0:round(seg_size/increased_p)-1
           seg_p_increased(i*increased_p+1) = 1;
       end
       for i = 1:length(AR_coeff)
           seg_p_increased = AR(seg_p_increased, AR_coeff(i, :));
       end
       
       % Synthesize the segment of voiced speech with lower pitch period
       decresed_p = round(pitch/pitch_modulation);
       seg_p_decresed = zeros(1, seg_size);
       for i = 0:round(seg_size/decresed_p)-1
           seg_p_decresed(i*decresed_p+1) = 1;
       end
       for i = 1:length(AR_coeff)
           seg_p_decresed = AR(seg_p_decresed, AR_coeff(i, :));
       end
    end
%   Make the segment energy equal to that of originally recorded segment
    seg = seg*sqrt(energy/(seg*seg'));
    seg_p_increased = seg_p_increased*sqrt(energy/(seg_p_increased*seg_p_increased'));
    seg_p_decresed = seg_p_decresed*sqrt(energy/(seg_p_decresed*seg_p_decresed'));
end

end

function y_ar = AR(x,a)
%   AR filter as done in class
    na = length(a);
    mem = zeros(1,na);
    for i = 1:length(x)
        % Shift memory and add new value
        y_ar(i) = a*mem' + x(i);
        mem(2:na) = mem(1:na-1);
        mem(1) = y_ar(i);
    end
end

%calculate the autocorrelation value up to maxdelta
function RX=Acor(X,maxDelta)
nx = length(X);
RX = zeros(1,maxDelta+1);
for idel = 0:maxDelta
    RX(idel+1) = X(1:nx-idel)'*X(idel+1:nx);
end
end

function signal = record()
recObj = audiorecorder(8000,8,1);   % sets mic as input
disp('Start speaking now')          % prompt speaker
recordblocking(recObj, 2);          % records from mic for 2 sec
disp('End of recording');           % indicate end
signal = getaudiodata(recObj)';     % write data in real-valued array
end

% Load dataset from file
function data = fromFile(filename)
fid = fopen(filename,'r');
data = fscanf(fid,'%f\n');
fclose(fid);
end

% Save a dataset to a file
function toFile(data, filename)
fid = fopen(filename,'w');
for i = 1:length(data)
    fprintf(fid,'%6.2f\n',data(i)); % write using "," separator
end
fclose(fid);
end
