%% HC-05 Bluetooth COM Port Scanner

% User Settings 
baudRate = 9600;
portRange = 3:20;		% COM3 to COM20
maxReadAttempts = 3;	% Number of read attempts per port

fprintf('\nScanning COM ports for HC-05 Bluetooth module...\n\n');

validPorts = [];

for i = portRange
    port = "COM" + i;
    fprintf('Checking %s... ', port);

    try
        bt = serialport(port, baudRate, "Timeout", 3);
        configureTerminator(bt, "LF");
        pause(2); % Allow port to stabilize

        isValid = false;
        for j = 1:maxReadAttempts
            try
                line = readline(bt);
                val = str2double(strtrim(line));
                if ~isnan(val) && val >= 0 && val <= 1023
                    fprintf('Valid ECG value received: %g\n', val);
                    validPorts(end+1) = i;
                    isValid = true;
                    break;
                end
            catch
                % Read error; skip attempt
            end
        end
        clear bt;

        if ~isValid
            fprintf('No valid data received.\n');
        end

    catch ME
        fprintf('Failed to open port (%s)\n', ME.message);
    end
end

% Summary
if isempty(validPorts)
    fprintf('\nNo valid HC-05 Bluetooth ports detected.\n');
else
    fprintf('\nDetected working HC-05 Bluetooth COM ports:\n');
    fprintf('   COM%d\n', validPorts);
    fprintf('\nUse one of these ports in your ECG acquisition script.\n');
end