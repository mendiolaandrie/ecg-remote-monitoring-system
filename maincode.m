% 1. User Settings
port = "COM15";         % Your HC-05 Bluetooth COM port
baudRate = 9600;		% HC-05 default baud rate
Fs = 320;               % Sampling frequency from Arduino
runTimeSeconds = 60;    % How long to record
windowSize = 1000;		% Live plot buffer
fprintf('--- ECG Real-Time Wireless Acquisition and Filtering ---\n');
fprintf('Sampling rate: %d Hz\n', Fs);
fprintf('Bluetooth Port: %s\n', port);
% 2. Connect to HC-05 Bluetooth
disp("Connecting to HC-05 on " + port + " ...");
bt = serialport(port, baudRate);
configureTerminator(bt, "LF");
bt.Timeout = 5;
flush(bt);
pause(2);               % Let Bluetooth settle
disp("Connected!");
% 3. Live Raw Data Plot During Acquisition
dataBuffer = zeros(1, windowSize);
figure('Name', 'Live ECG Acquisition');
hPlot = plot(1:windowSize, dataBuffer, 'b');
ylim([0 1023]);
xlabel('Sample Index');
ylabel('ADC Value');
title('Live ECG: Raw Signal (Bluetooth)');
grid on;
disp("Starting real-time acquisition...");
startTime = tic;
allRawData = [];
while toc(startTime) < runTimeSeconds
   try
       val = str2double(strtrim(readline(bt)));
       if isnan(val) || val < 0 || val > 1023
           val = 0;
       end
   catch
       val = 0;
   end
   allRawData = [allRawData, val];
   dataBuffer = [dataBuffer(2:end), val];
   set(hPlot, 'YData', dataBuffer);
   drawnow limitrate;
end
clear bt;
disp("Acquisition finished. Bluetooth port closed.");
% 4. Rescale ADC Values to mV
ecg_mv = (allRawData / 1023) * 3.3;
ecg_mv = ecg_mv(:); % Ensure column
t = (0:length(ecg_mv)-1) / Fs;
figure;
plot(t, ecg_mv);
title('Raw Acquired ECG Signal');
xlabel('Time (s)');
ylabel('Voltage (approx mV)');
grid on;
% 5. Baseline Wander Removal (Wavelet Decomposition)
[C, L] = wavedec(ecg_mv, 9, 'bior3.7');
% Sum all detail components
details = cell2mat(cellfun(@(lvl) wrcoef('d', C, L, 'bior3.7', lvl), num2cell(1:9), 'UniformOutput', false));
baseline_removed = sum(details, 2);
figure;
subplot(2,1,1); plot(t, ecg_mv); title('Raw ECG Signal'); grid on;
subplot(2,1,2); plot(t, baseline_removed); title('After Baseline Wander Removal'); grid on;
% 6. 50 Hz Notch Filter
Fnotch = 50;
BW = 2;
[b_notch,a_notch] = iirnotch(Fnotch/(Fs/2), BW/(Fs/2));
notched = filter(b_notch, a_notch, baseline_removed);
figure;
plot(t, notched);
title('After 50Hz Notch Filter');
grid on;
% 7. Low-pass Filter to Remove High-Frequency Noise
d = fdesign.lowpass('Fp,Fst,Ap,Ast', 0.4, 0.5, 1, 80);
Hd = design(d,'equiripple');
lowpassed = filter(Hd, notched);
% Ensure numeric vector
lowpassed = lowpassed(:);
lowpassed(isnan(lowpassed)) = 0;
figure;
plot(t, lowpassed);
title('After Low-Pass Filtering');
grid on;
% 8. Wavelet Denoising (MODWT)
if length(lowpassed) < 2
   error("Signal too short for MODWT.");
end
wt = modwt(lowpassed, 4, 'sym4');
wtrec = zeros(size(wt));
wtrec(3:4,:) = wt(3:4,:);
denoised = imodwt(wtrec, 'sym4');
figure;
plot(t, denoised);
title('Final Denoised ECG Signal');
grid on;
% 9. PQRST Feature Detection
[Rpeaks,locs_r] = findpeaks(denoised, 'MinPeakHeight', 0.2*max(denoised), 'MinPeakDistance', round(0.6*Fs));
nohb_r = length(locs_r);
% S Peaks (local minima after R)
Speaks = zeros(size(locs_r));
locs_s = zeros(size(locs_r));
for i = 1:nohb_r
   range = locs_r(i):min(length(denoised), locs_r(i)+15);
   [Speaks(i),idx] = min(denoised(range));
   locs_s(i) = locs_r(i)+idx-1;
end
% Q Peaks (local minima before R)
Qpeaks = zeros(size(locs_r));
locs_q = zeros(size(locs_r));
for i = 1:nohb_r
   range = max(1, locs_r(i)-15):locs_r(i);
   [Qpeaks(i),idx] = min(denoised(range));
   locs_q(i) = range(1)+idx-1;
end
% P Peaks (before Q)
Ppeaks = zeros(size(locs_r));
locs_p = zeros(size(locs_r));
for i = 1:nohb_r
   range = max(1, locs_q(i)-60):locs_q(i);
   [Ppeaks(i),idx] = max(denoised(range));
   locs_p(i) = range(1)+idx-1;
end
% T Peaks (after S)
Tpeaks = [];
locs_t = [];
for i = 1:nohb_r
   if locs_s(i)+130 <= length(denoised)
       range = locs_s(i):locs_s(i)+130;
       [tp,tp_idx] = max(denoised(range));
       Tpeaks(end+1) = tp;
       locs_t(end+1) = locs_s(i)+tp_idx-1;
   end
end
% 10. Heart Rate Calculation
timelimit = length(denoised)/Fs;
hbpermin = (nohb_r * 60) / timelimit;
fprintf('\nEstimated Heart Rate: %.1f BPM\n', hbpermin);
% 11. Plot with PQRST Markers
figure;
plot(t, denoised);
grid on;
hold on;
plot(t(locs_r), Rpeaks, '^r', 'DisplayName', 'R peaks');
plot(t(locs_s), Speaks, 'xm', 'DisplayName', 'S peaks');
plot(t(locs_q), Qpeaks, '*g', 'DisplayName', 'Q peaks');
plot(t(locs_p), Ppeaks, 'om', 'DisplayName', 'P peaks');
plot(t(locs_t), Tpeaks, 'pk', 'DisplayName', 'T peaks');
xlabel('Time (s)');
ylabel('Voltage (mV approx)');
title(sprintf('PQRST Detection and Heart Rate: %.1f BPM', hbpermin));
legend;
% 12. Clinical Feature Extraction and Summary
RR_intervals_sec = diff(locs_r) / Fs;
mean_RR_ms = mean(RR_intervals_sec) * 1000;
PR_intervals_sec = (locs_r - locs_p) / Fs;
mean_PR_ms = mean(PR_intervals_sec) * 1000;
QRS_durations_sec = (locs_s - locs_q) / Fs;
mean_QRS_ms = mean(QRS_durations_sec) * 1000;
% QT intervals
numQT = min(length(locs_r), length(locs_t));
QT_intervals_sec = (locs_t(1:numQT) - locs_r(1:numQT)) / Fs;
mean_QT_ms = mean(QT_intervals_sec) * 1000;
% QTc (Bazett correction)
QTc_ms = QT_intervals_sec ./ sqrt(RR_intervals_sec(1:numQT));
mean_QTc_ms = mean(QTc_ms)*1000;
fprintf('\n--- Clinical Feature Summary ---\n');
fprintf('Mean Heart Rate: %.1f BPM\n', hbpermin);
fprintf('Mean RR Interval: %.1f ms\n', mean_RR_ms);
fprintf('Mean PR Interval: %.1f ms\n', mean_PR_ms);
fprintf('Mean QRS Duration: %.1f ms\n', mean_QRS_ms);
fprintf('Mean QT Interval: %.1f ms\n', mean_QT_ms);
fprintf('Mean QTc Interval: %.1f ms\n', mean_QTc_ms);
fprintf('-------------------------------\n');
