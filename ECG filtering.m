close all
y = load('noisyECG.mat');
n = length(y.y);
fs = 250;
T = 1/fs;
L = length(y.y);
w0 = (0:n-1);
figure
plot(y.y); 
title ('ECG signal with noise');
xlabel('t(1/250s)');
ylabel('ECG signal');





s = fftshift(fft(y.y));
w1 = linspace(-fs/2,fs/2-fs/n,n);
figure
plot(w1,abs(s));
title ('fft and fftshift response');
xlabel('frequency(Hz)');
ylabel('Magnitude');


%second order IIR notch filter
fn = 60;
wn = 2*pi*fn/fs;
r = 0.95;
b = [1 -2*cos(wn) 1];
a = [1 -2*r*cos(wn) r.^2];


[H, W] = freqz(b,a);
b = b/abs(H(1)); 
figure
subplot(2,1,1);
plot(y.y);
title ('Original data');
xlabel('t(1/250s)');
ylabel('(ECGsignal)');

subplot(2,1,2);
plot(firstResult);
title ('After IIR notch filtering');
xlabel('t(1/250s)');
ylabel('(ECGsignal)');

figure
subplot(2,1,1);
plot(w1,abs(s));
title ('Original fft and fftshift response');
xlabel('frequency(Hz)');
ylabel('Magnitude');

subplot(2,1,2);
plot(w1,abs(fftshift(fft(firstResult))));
title ('After IIR notch filtering');
xlabel('frequency(Hz)');
ylabel('Magnitude');


%Chebyshev II lowpass filter
wp = 42/fs;
ws = 60/fs;
Rp = 1;
As = 50;

[N1, W1] = cheb2ord(wp, ws, Rp, As);


[b1, a1] = cheby2(N1, As, W1);
[H1, W1]=freqz(b1,a1);
fprintf('The order of the transfer function of the filter: %d\n',N1);

fprintf('The top coefficients of the transfer function of the filter: ');
for i = 1:length(b1)
    fprintf('%.5g ',b1(i));
end
fprintf('\n');
fprintf('The bottom coefficients of the transfer function of the filter: ');
for i = 1:length(a1)
    fprintf('%.5g ',a1(i));
end
fprintf('\n');
% plot(abs(fft(firstResult)))
secondResult = filter(b1,a1,firstResult);

figure
subplot(2,1,1)
plot(W1/pi, 20*log10(abs(H1)));
title ('Chebyshev II lowpass filter');
xlabel('frequency(Hz)');
ylabel('Magnitude of the frequency response');
subplot(2,1,2)
zplane(b1,a1);

figure
subplot(2,1,1)

plot(firstResult);
title ('After IIR notch filtering');
xlabel('t(1/250s)');
ylabel('(ECGsignal)');
subplot(2,1,2)
plot(secondResult);
title ('After Chebyshev II lowpass filter ');
xlabel('t(1/250s)');
ylabel('(ECGsignal)');


%Chebyshev II high pass filter
wp=0.8;
ws=0.75;
Rp=1;
As=50;
[N2, W2] = cheb2ord(wp, ws, Rp, As);
[b2, a2] = cheby2(N2, As, W2,'high');

[H2, W2]=freqz(b2,a2);
% plot(W2/pi, 20*log10(abs(H2)));
fprintf('The order of the transfer function of the filter: %d\n',N2);

fprintf('The top coefficients of the transfer function of the filter: ');
for i = 1:length(b2) 
    fprintf('%.5g ',b2(i));
end
fprintf('\n');
fprintf('The bottom coefficients of the transfer function of the filter: ');
for i = 1:length(a2)
    fprintf('%.5g ',a2(i));
end
fprintf('\n');
thirdResult = filter(b2,a2,secondResult);

figure
subplot(2,1,1)
plot(secondResult);
title ('After Chebyshev II lowpass filter');
xlabel('t(1/250s)');
ylabel('(ECGsignal)');
subplot(2,1,2)
plot(thirdResult);
title ('After Chebyshev II high pass filter ');
xlabel('t(1/250s)');
ylabel('(ECGsignal)');


maxValue = max(thirdResult);
threshold = 0.7*maxValue;

%Thresholding filtered result
for i = 1:n
    if(thirdResult(i)<threshold)
        thirdResult(i) = 0;
    else
        thirdResult(i) = thirdResult(i) - threshold;
    end
end
[pks,locs] = findpeaks(thirdResult);
figure
findpeaks(thirdResult);

for i = 1:length(pks)
    for j = i:length(pks)
        if((locs(j)-locs(i))<50&&i~=j)      
            if(pks(i) < pks(j))
                pks(i) = 0;
            else
                pks(j) = 0;
            end
            
        end
    end
    
end

count = 0;
var = 0;
for i = 2:length(pks)
    for j = i:length(pks)
        if(pks(i)~=0&&i~=j)
            var = var + abs(pks(i)-pks(i-1));
            count = count +1;
            break
        end
    end

end
avgVar = var/count;

heartRate = 0;
count = 0;
for i = 1:length(pks)
    for j = i:length(pks)
        if((locs(j)-locs(i))>50&&i~=j)      
            heartRate = heartRate + 60/((locs(j)-locs(i))/250);
            count = count +1;
            break
        end
    end
end
avgHeartRate = heartRate/count;
% 
% 
% count = 0;
% maxCurr = 0;
% rate = 0;
% countC = 0;
% var = 0;
% for i = firstValidIndex:n
%     
%     if(thirdResult(i)>0&& count>=50&&countC>0)
%         tmpR = 60/(count/250);
%         rate = rate + tmpR;
%         count = 0;
%         countC = countC +1;
%         tmpVar = abs(thirdResult(i)-maxCurr);
%         var = var+tmpVar;
%        
%         
%     else
%         if(thirdResult(i)>0&& count>=50&&countC==0)
%             countC = countC +1;
%         end
%     end
%     
%     
%     count = count + 1;
% end
format shortG
fprintf('average heart beat is: %.5g\n',avgHeartRate);
fprintf('average heart beat variability is: %.5g\n',avgVar);


