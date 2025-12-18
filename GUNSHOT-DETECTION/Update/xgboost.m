% 1. Setup Data Paths
% UPDATE THIS PATH to your actual folder location containing 'Gunshot' and 'NonGunshot'
dataFolder = 'Data'; 

% Create Datastore
ads = audioDatastore(dataFolder, 'IncludeSubfolders', true, 'LabelSource', 'foldernames');
[trainDatastore, testDatastore] = splitEachLabel(ads, 0.8);

%% 2. Define Feature Extraction (Synced with Simulink)
fs = 44100;
windowLength = 1024;
overlapLength = 512;

afe = audioFeatureExtractor('SampleRate', fs, ...
    'Window', hann(windowLength, 'periodic'), ...
    'OverlapLength', overlapLength, ...
    'mfcc', true);

% ENABLE Delta to match your Simulink "Append delta" checkbox
afe.mfccDelta = true; 
% DISABLE DeltaDelta (Assuming you only checked "Append delta")
afe.mfccDeltaDelta = false;

%% 3. Extract Features (Training Set)
disp('Extracting features from Training Set...');
trainFeatures = [];
trainLabels = [];

reset(trainDatastore);
while hasdata(trainDatastore)
    [audioIn, info] = read(trainDatastore);
    if size(audioIn, 2) > 1, audioIn = mean(audioIn, 2); end % Mono
    
    features = extract(afe, audioIn);
    trainFeatures = [trainFeatures; features];
    
    numFrames = size(features, 1);
    trainLabels = [trainLabels; repmat(info.Label, numFrames, 1)];
end

%% 4. Train Model
disp('Training XGBoost/Ensemble Model...');
trainedModel = fitcensemble(trainFeatures, trainLabels, ...
    'Method', 'Bag', ...
    'NumLearningCycles', 50, ...
    'Learners', templateTree('MaxNumSplits', 20));

%% 5. Calculate Accuracy (Test Set)
disp('Testing Model Accuracy...');
testFeatures = [];
testLabels = [];

reset(testDatastore);
while hasdata(testDatastore)
    [audioIn, info] = read(testDatastore);
    if size(audioIn, 2) > 1, audioIn = mean(audioIn, 2); end
    
    features = extract(afe, audioIn);
    testFeatures = [testFeatures; features];
    
    numFrames = size(features, 1);
    testLabels = [testLabels; repmat(info.Label, numFrames, 1)];
end

% Predict on test data
[predictedLabels, scores] = predict(trainedModel, testFeatures);

% Calculate Percentage
accuracy = sum(predictedLabels == testLabels) / length(testLabels) * 100;
fprintf('<strong>Final Model Accuracy: %.2f%%</strong>\n', accuracy);

% Plot Confusion Matrix
figure;
confusionchart(testLabels, predictedLabels);
title(['Confusion Matrix (Accuracy: ' num2str(accuracy, '%.1f') '%)']);

%% 6. Save for Simulink
save('gunshotModel.mat', 'trainedModel');
disp('Model saved! Ready for Simulink.');