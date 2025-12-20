%% 1. Setup & Initial Scan
dataFolder = 'Data'; % Your folder with 'Gunshot' and 'NonGunshot'
ads = audioDatastore(dataFolder, 'IncludeSubfolders', true, 'LabelSource', 'foldernames');

disp('--- Initial Dataset Status ---');
labelCounts = countEachLabel(ads);
disp(labelCounts);

% Find the target count (Match the largest class)
maxCount = max(labelCounts.Count);
uniqueLabels = labelCounts.Label;

%% 2. STEP A: BALANCE THE DATASET (Make Equal)
disp('--- Step A: Balancing Dataset (Filling gaps) ---');

for i = 1:length(uniqueLabels)
    thisLabel = uniqueLabels(i);
    thisCount = labelCounts.Count(labelCounts.Label == thisLabel);
    
    if thisCount < maxCount
        numNeeded = maxCount - thisCount;
        fprintf('Class "%s" needs %d files to balance. Generating...\n', string(thisLabel), numNeeded);
        generateAugmentedFiles(ads, thisLabel, numNeeded, '_balance');
    else
        fprintf('Class "%s" is already the max size. Skipping balance.\n', string(thisLabel));
    end
end

%% 3. STEP B: DOUBLE THE DATASET (Expand All)
% Now that they are equal, we force-add more data to BOTH classes
disp('--- Step B: Doubling the Entire Dataset ---');

% Re-scan to get the balanced list
ads = audioDatastore(dataFolder, 'IncludeSubfolders', true, 'LabelSource', 'foldernames');
labelCounts = countEachLabel(ads);

for i = 1:length(uniqueLabels)
    thisLabel = uniqueLabels(i);
    currentCount = labelCounts.Count(labelCounts.Label == thisLabel);
    
    % We want to add an amount equal to the current count (Doubling it)
    fprintf('Doubling class "%s" (Adding %d more files)...\n', string(thisLabel), currentCount);
    generateAugmentedFiles(ads, thisLabel, currentCount, '_double');
end

disp('--- Data Prep Complete. Re-scanning for Training... ---');



%% %% --- HELPER 2: File Generator (Corrected to fix Clipping) ---
function generateAugmentedFiles(adsFull, targetLabel, numNeeded, suffix)
    % Filter to only the files of the target label
    subAds = subset(adsFull, adsFull.Labels == targetLabel);
    files = subAds.Files;
    
    generated = 0;
    fileIdx = 1;
    
    while generated < numNeeded
        % Cycle through original files
        originalFile = files{mod(fileIdx-1, length(files)) + 1};
        [audioIn, fs] = audioread(originalFile);
        
        % Determine augmentation type
        method = mod(generated, 3);
        
        if method == 0 % White Noise
            noise = 0.015 * randn(size(audioIn)); 
            audioAug = audioIn + noise;
            augType = '_noise';
            
        elseif method == 1 % Pitch Shift
            semitones = (rand * 4) - 2; 
            audioAug = shiftPitch(audioIn, semitones);
            augType = '_pitch';
            
        else % Time Shift
            shiftAmount = round(0.2 * fs); 
            audioAug = circshift(audioIn, shiftAmount);
            augType = '_time';
        end
        
        % --- FIX FOR CLIPPING (Normalize) ---
        maxVal = max(abs(audioAug));
        if maxVal > 1
            audioAug = audioAug / maxVal; % Scale down so max is exactly 1.0
        end
        % ------------------------------------
        
        % Save
        [folderPath, name, ext] = fileparts(originalFile);
        newFileName = fullfile(folderPath, [name suffix '_' num2str(generated) augType ext]);
        
        if ~isfile(newFileName)
            audiowrite(newFileName, audioAug, fs);
            generated = generated + 1;
        end
        
        fileIdx = fileIdx + 1;
    end
end
