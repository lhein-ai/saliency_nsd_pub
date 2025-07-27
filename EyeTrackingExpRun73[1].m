
%% Receiving ground-truth data of visual saliency via Eye-Tracking of pictures of the Microsoft Coco Images Dataset
%% that were presented to the participants in Natrual Senses Dataset (NSD-Stimuli)

% This code runs an EyeTracking Experiment with EyeLink 1000 and shows
% images that were presented to the participants in the Natural Senses
% Dataset (N = 8) by Allen et al (2021). These images were retrieved by the
% Microsoft Coco images Dataset and show random pictures of daily life
% (e.g. landscapes, people, animals, drinks & food, sports). The shared1000
% present images that were presented to all NSD participants and are shown
% to all patricipants within this experiment as well. (Moreover each
% particpant gets presented 100 pictures of the other 72.000 images that were
% shown to one of the particpants with in the NSD experiment.)
% To get the particpants attention and gain for higher comparison of the
% data retrived from Allen et al. (2021), the subjects are asked to perform
% a continous recognition task detecting old_images (images that have been
% preseted to them before) correctly.

% The probability that an image that is shown is an old one is set to 1/4.
% The repetition around the blocks is equalized so that the particpants
% can't just press new_image all the time or guess and retrieve a high
% feedback score.
% The feedback to the participants gets presented after every block and
% reflects their ability to correctly detect the old images as 'old' and the
% 'new' images as 'new'.

%% Function languagetransfer: gets the desired text language (English or German)
% and presents the text parts in the desired language. This function is not
% a local function, so it is ectern.

%% Local function: get_img_seq
% gets the image_seq with desired repetition probability distrubted the
% same during the whole experiment.

%%% Written by Lisa Heinemann
%%% Last edited: 03/01/2024

%% TODO 1: Run function
function EyeTrackingExpRun73

% set screen resolution to 1920 x 1080 to control experimental screen
win_horizontal = 1920;
win_vertical = 1080;

%----------------------------------------------------------------------
%                     Open experimental screen
%----------------------------------------------------------------------
%% TODO 3: Change to 0 so that they get edited
Screen('Preference', 'SkipSyncTests', 0);

% Call of some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);

% ----------------------------------------------------------------------
%                           Base directory
%----------------------------------------------------------------------

%Direction to Apple (MacOS) - Windows - Linux (Fedora) - Switch Case
%changes directions automatically

switch computer()
    %% Linux
    case 'GLNXA64'
        base_dir = '/home/lab/Documents/lheinemann-ET-NSD-cocoImg';

        %% Apple (macOS)
    case 'MACI64'
        %base_dir = '/Volumes/bartels_data/lheinemann';

        base_dir = '/Users/lisaheinemann/Documents/MATLAB';

        %% Windows
    case 'PCWIN64'
        base_dir= 'C:/Users/Lisa/MATLAB Drive/Matlab';
end
%% TODO 4: Use of eye link, yes or no? Press 1 to activate 0 means not connected
% Needs to be zero so that code runs on computer without eye tracking 
use_of_eyelink = 1;

%_________________________________________________________________________
% Load pictures for test run

% cell array filnamesC for loading stimuli, 7 pictures, 3 of them are
% repeated, pictures retrieved from Wikimedia, license-free usage and
% cropped seantically and rescaled to same resolution as NSD pic (425x425
% pixels)
filenamesC = {'CocoAlikeImg1c.jpg';
    'CocoAlikeImg2c.jpg';
    'CocoAlikeImg1c.jpg';
    'CocoAlikeImg3c.jpg';
    'CocoAlikeImg4c.jpg';
    'CocoAlikeImg5c.jpg';
    'CocoAlikeImg6c.jpg';
    'CocoAlikeImg3c.jpg';
    'CocoAlikeImg7c.jpg';
    'CocoAlikeImg4c.jpg';
    };
%----------------------------------------------------------------------
%               Constant params used in the Experiment
%----------------------------------------------------------------------
%% TODO 5: Change to the desired amount of images

%% number of different images to present in main experiment (all blocks)
% E.g. with n_imgs = 16, n_blocks = 4 and old_prob = 1/3 this makes 24
% pictures presented in the whole experiment and 6 per trial

%n_imgs = 1000
% to test: n_imgs = 24, n_blocks = 4, old_prob = 1/4
n_imgs = 1000;

%1050 or 999 best number of pictures because 1050/(1/4)-1) =1400 or 960

%% How many blocks do we want in total? (usually 5?) --> 1050/5 = 210
n_blocks = 5;

% What fraction of trials do we want to require an "old" response?
old_prob = 1/5;

%%%%%%%%%%%%%%%%%%%%%%%%%% timing parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Duratin pictures are shown
requested_stimulusDuration = 3; % in seconds

% Duration of Interstimulus interval (ISI) - blank space with fixation dot
requested_fixationDuration = 1;  % xxx interval time in seconds

%%%%%%%%%%%%%%%%%% Initalize the data structure %%%%%%%%%%%%%%%%%%%%%%%%
% Data saves the whole information of experiment: details about the
% participant, Screen & Timing Information (Stim-on and Stim-offset) as
% well as information of pictures, there should be 4 files per trial
data = struct();

% In "testRunImages" the information relevant for each picture of test run gets saved
data.testRunImages = struct();

% In "pictures" the information relevant for each picture gets saved
data.images = struct();
%----------------------------------------------------------------------
%    Open PTB-3 to prepare screen and get WindowSize (pixel dimensions)
%----------------------------------------------------------------------
%%%%%%%%%%%%%%%%% screenfunction and preferences %%%%%%%%%%%%%%%%%%%%%%
% Get the screen number
screenVector = Screen('Screens');
%% Adjust screenNumbers automatically
% If screenvector = one screen then PTB3-experiment is presenetd on it with transparent
% window, if screenVector (0,1) or (0,1,2) experiment runs on the screen
% with the highst number minus one: so in the case of (0,1,2) on 1 = screen
% 2

if (max(screenVector > 0))
    % screenNumber = max(screenVector) - 1 to reach the second monitor
    screenNumber = max(screenVector) - 1;
else
    % if just one screen is connected, screenVector = 0
    screenNumber = max(screenVector);
end

% If problems exist debugging ideas:

%screenNumber = min(screenVector)+1;
%screenNumber = 1;

% Define black, white and grey
black = BlackIndex(screenNumber);
white = WhiteIndex(screenNumber);
grey = white / 2;

% Prepare imaging pipeline for some fitting work:
PsychImaging('PrepareConfiguration');

% Set display resolution
PsychImaging('AddTask', 'General', 'UsePanelFitter', ...
    [win_horizontal, win_vertical], 'Aspect');

%% Make the screen transparent if just one screen is used
PsychDebugWindowConfiguration;

% Open an on screen window and color it grey
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, grey);

% Flip to clear
Screen('Flip', window);

%S et blend function for the screen
Screen('BlendFunction', window, 'GL_SRC_ALPHA','GL_ONE_MINUS_SRC_ALPHA');

% Query the frame duration (inter frame interval)
ifi = Screen('GetFlipInterval', window);

% Set the text size
Screen('TextSize', window, 20);

% Get the size of the on screen window in pixels
[screenXpixels, screenYpixels] = Screen('WindowSize', window);

% Get the centre coordinate of the window in pixels
% For help see: help RectCenter
[xCenter, yCenter] = RectCenter(windowRect);

%%%%%%%%%%%%%%%%%%%%%%%% Screen details %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TODO: Change if needed

% Screen distance (in centimeters) in exp 55

%% TODO 6: Adapt to experiment
% screen dimensions (in centimeters)
data.screen_dist = 55;

%% TODO 7: Suit to experimental monitor (width 60 cm and height 34 cm)
% Used screen in this experiment
% Eizo Flex Scan EV2736W

data.screen_width = 60;
data.screen_height = 34;

% Calculate screen diagonal (in cm) and  write it to log file
data.screen_diag = sqrt(data.screen_width^2+ data.screen_height^2);

%%%%%%%%%%%%%%% screen resolution (in centimeters) %%%%%%%%%%%%%%%%%%%%%

% information of pixel size of monitor from
% https://www.displayspecifications.com/de/model/7ea254
% Code gets the amount of pixels directly from the Experimental Screen

% Amount of pixels on the vertical line
data.screen_pix_hor = screenXpixels;

% Amount of pixels on the horizontal line
data.screen_pix_ver = screenYpixels;

%%%%%%%%%%%%%%%%%%%%%%%%% Calculate PPCM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate PPCM (Pixels per centimeter)
data.ppcm = sqrt(data.screen_pix_hor^2 + data.screen_pix_ver^2) / data.screen_diag;

% Display the result and write the Pixels per centimeter to the log file
fprintf('Pixels per cm of screen (ppcm) = %.2f.\n', data.ppcm);

% Prompt for the visual angle in degrees
data.visualAngle = 8.4;

%%%%%%%%%%%%%%% Calculate Pixel size given the angle %%%%%%%%%%%%%%%%%%%%

% Calculate pixel size for the given visual angle
data.cm_sz = 2 * data.screen_dist * tan(deg2rad(data.visualAngle) / 2);

data.pix_siz = 2 * data.screen_dist * tan(deg2rad(data.visualAngle) / 2) * data.ppcm;

% Display the calculated pixel size and write the calculated adapted pixel size to the log file
fprintf('Pixel size for given visual angle: %.2f pixels\n', data.pix_siz);
fprintf('Adapted Pixel size: %.1f\n', data.pix_siz);

%%%%%%%%%%%%%%%%%%%%%%%%% Subject data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Welcome to our experiment');

prompt = "Subject number: ";
data.nr = input(prompt);

prompt = "Age? ";
data.age = input(prompt);

% Controlling that input of gender is only 'm', 'f', or 'd'
prompt = "Gender? m/f/d ";
while true
    data.gender = input(prompt, 's');

    if ismember(data.gender, {'m', 'f', 'd'})
        break
    else
        prompt = "Please enter only 'm', 'f', or 'd': ";
    end
end

% Prompt for block number
prompt = "current block number? Enter 1,2,3,4 or 5!";
data.currentBlockNumber = input(prompt);

% Prompt for language
prompt = "current language? English or German!";
data.currentLanguage = input(prompt,"s");

% get the date for file saving
todays_date = yyyymmdd(datetime());
data.date = todays_date;

% function call language transfer to switch languages
[text] = languagetransfer(data.currentLanguage, data.currentBlockNumber, n_blocks);

%%%%%%%%%%%%%%%%%%%%%%%% EYELINK Part I %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create edf file of eyelink
if use_of_eyelink
    edfFile = ['ET_', num2str(data.nr), '_', num2str(data.currentBlockNumber)];
end

%%%%%%%%%%%%%%%%%%%%%% for randomization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%reset the random number generator so that at each trial depending on
%subject no. different pictures are shown
rng(data.nr);

%----------------------------------------------------------------------
%                     Load stimuli
%----------------------------------------------------------------------

% Get path of stimuli in base driectory
img_dir = fullfile(base_dir, 'shared1000');
addpath(img_dir);

% Read stimuli
pictures = dir(fullfile(img_dir, '*.jpg'));

% Column showing pictures name
pictures = transpose({pictures.name});


%%%%%%%%%%%%%%%%%%%%%%%% EYELINK Part II %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------------------------------------------------------
%                     Eye-Tracking Setup
%----------------------------------------------------------------------
% Open eye tracker on screen to record data
if use_of_eyelink
    EyelinkInit(0,1);
    el = EyelinkInitDefaults(window);

    [~, vs] = Eyelink('GetTrackerVersion');
    fprintf('Runnning experiment on a ''%s'' tracker.\n', vs);

    % factor of reduced display and thus calibration area 
    Eyelink('command', 'calibration_area_proportion = 0.3 0.3');
    Eyelink('command', 'validation_area_proportion = 0.3 0.3');

    Eyelink('command', 'file_event_filter = LEFT, RIGHT, FIXATION, SACCADE, BLINK, MESSAGE');
    Eyelink('command', 'file_sample_data = LEFT,RIGHT,GAZE,AREA');

    %set link data (used for gaze cursor)
    Eyelink('command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
    % make sure that we get gaze data from the Eyelink
    Eyelink('Command', 'link_sample_data = LEFT,RIGHT,GAZE,AREA');

    Eyelink('command', 'screen_pixel_coords = %ld %ld %ld %ld', 0, 0, screenXpixels-1, screenYpixels-1);
    Eyelink('message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, screenXpixels-1, screenYpixels-1);
    % set calibration type
    Eyelink('command', 'calibration_type = HV13');

    %%% FIXME Uncomment if needed
    %proof remaining connection of Eyelink
    if Eyelink('IsConnected')~= 1
        fprintf('not connected, clean up\n');
        Screen('CloseAll');
        return;
    end
end
%----------------------------------------------------------------------
%               Start of experiment with Psychtoolbox
%----------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%% dependent on block number %%%%%%%%%%%%%%%%%%%%%%%

try

%%%%%%%%%%%%%%%%%%%%%%%% EYELINK Part III %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if eyetracker is used in experiment
if use_of_eyelink
    % status if edf (= eyetracking file) is opened
    % returns 0 if success, else error code
    status = Eyelink('Openfile', edfFile);
    if status ~=0
        fprintf('Cannot create EDF file ''%s'' ', edfFile);
        %Eyelink('Shutdown');
        Screen('CloseAll');
        return;
    end
end

%%%%%%%%%%% Start with instructions of the experiment %%%%%%%%%%%%%%%

if data.currentBlockNumber == 1

    DrawFormattedText(window, text.welcomeText, 'center', 'center', white);

    % Flip to the screen
    Screen('Flip', window);

    % Now we have drawn to the screen we wait for a keyboard button press (any
    % key) to terminate the demo
    KbStrokeWait;

    % Switches to next Page ScreenDuringExperimenth

    %%  Draw all the text in one go
    DrawFormattedText(window, text.introductionText1 ,'center', 'center', white);

    % Flip to the screen
    Screen('Flip', window);

    % Now we have drawn to the screen we wait for a keyboard button press (any
    % key) to terminate the demo
    KbStrokeWait;

    % Switches to next page: screen during testrun of experiment instructions


    %% Draw all the text in one go
    DrawFormattedText(window, text.introductionText2 ,'center', 'center', white);

    % Flip to the screen
    Screen('Flip', window);

    % Now we have drawn to the screen we wait for a keyboard button press (any
    % key) to terminate the demo
    KbStrokeWait;
    % %------------------------------------------------------------------
    % %                           Start with test run
    % %------------------------------------------------------------------

    % show first two pictures of NSD folder (nsd_stimulus_id_2951 and
    % nsd_stimulus_id_2991) show first one twice to see whether participants
    % understand the task and press 'j' while presenting third picture for
    % 'old picture'

    % Access the content of file enamesC-Array using indexing
    element1 = filenamesC{1};  % Access the 1st element, 1st pic - nsd_stimulus_id_2951.jpg
    element2 = filenamesC{2};  % Access the 2nd element, 2nd pic - nsd_stimulus_id_2991.jpg
    element3 = filenamesC{3};  % Access the 3rd element, 1st pic again --> repetition 1: nsd_stimulus_id_2951.jpg
    element4 = filenamesC{4};  % Access the 4th element, 3rd pic - nsd_stimulus_id_3050.jpg
    element5 = filenamesC{5};  % Access the 5th element, 4th pic - nsd_stimulus_id_3078.jpg
    element6 = filenamesC{6};  % Access the 6th element, 5th pic - nsd_stimulus_id_3147.jpg
    element7 = filenamesC{7};  % Access the 7th element, 6th pic - nsd_stimulus_id_3158.jpg
    element8 = filenamesC{8};  % Access the 8th element, 3rd pic again --> rep.: 2 nsd_stimulus_id_3050.jpg
    element9 = filenamesC{9};  % Access the 9th element, 7th pic- nsd_stimulus_id_3165.jpg
    element10 = filenamesC{10};  % Access the 10th element, 4th pic again --> rep.:3 nsd_stimulus_id_3078.jpg


    % Display the content
    disp(element1);
    disp(element2);
    disp(element3);
    disp(element4);
    disp(element5);
    disp(element6);
    disp(element7);
    disp(element8);
    disp(element9);
    disp(element10);


    % definition target size of resized image to fit given visual angle
    % of 8.4
    target_size = [data.pix_siz, data.pix_siz];

    n_practice_trials = 10;
    trial_dur = requested_fixationDuration + requested_stimulusDuration;

    req_img_onsets_TR = (0:(n_practice_trials-1)) * trial_dur + requested_fixationDuration;
    req_img_offsets_TR = req_img_onsets_TR + requested_stimulusDuration;

    %-----------------------------------------------------------------------
    %            TR - ISI: Draw fixition dot in blank space (1s) - in
    %            prpartion prior to loop
    %-----------------------------------------------------------------------

    % prepare drawing dots
    testrun_start_time = Screen('Flip', window);
    Screen('DrawDots', window, [xCenter; yCenter], 10, black, [], 2);

    % Save trial information in advance to data log
    data.block_start_time_testr = testrun_start_time;
    data.req_img_onsets_testr = req_img_onsets_TR;
    data.req_img_offsets_testr= req_img_offsets_TR;

    for j = 1:n_practice_trials

        %Get the presentation time of the ISI
        data.testRunImages(j).isi_on = Screen('Flip', window, requested_fixationDuration)-testrun_start_time; %starttime ISI

        data.testRunImages(j).isi_off = WaitSecs(1)-testrun_start_time; % endtime ISI

        % isRecDuration = complete time interval ISI was shown (supposed to be 1s)
        data.testRunImages(j).isiRecDuration = data.testRunImages(j).isi_off - data.testRunImages(j).isi_on;

        %-----------------------------------------------------------------------
        %                Load stimuli  - TR
        %-----------------------------------------------------------------------

        % Load stimuli (images)
        % Load stimulus image based on filenames from filenamesC
        currentImageFilename = filenamesC{j};
        currentImage = imread(currentImageFilename);

        % Resize the image
        resized_image = imresize(currentImage, target_size);

        % Present the resized stimulus
        texture = Screen('MakeTexture', window, resized_image);

        %------------------------------------------------------------------
        %            Test run: Draw Images
        %------------------------------------------------------------------
        % Draw the texture and flip the window
        Screen('DrawTexture', window, texture);
        % data.images(j).T_stim_on_t  = Screen('Flip', window, data.images(j).req_img_onset); % add req_img_onset;
        img_onset_t = Screen('Flip', window, req_img_onsets_TR(j) + testrun_start_time - .5 * ifi);

        % T_t1 = 0; % Reset the time counter
        old_responses = []; % Iniialize to empty vector to be filled with key presses
        RTs = []; % Initialize to empty vector with response times relative to image onset


        %% detect key presses during while loop and see if it's time to break the loop
        response_registered = false;
        correct_response = false;
        hits = false;
        misses = false;
        correctRejections = false; 
        falseAlarms = false; 

        while GetSecs < testrun_start_time + req_img_offsets_TR(j) - 1.5 * ifi

            [keyIsDown, secs, keyCode, ~] = KbCheck;
            if keyIsDown
                % abort experiment in pressing x during execution of while loop
                if keyCode(KbName('x'))
                    ShowCursor;
                    sca;
                    fclose all;
                    Screen('CloseAll');
                    return
                end
                if (keyCode(KbName('4')) || keyCode(KbName('5'))) && response_registered==false
                    RTs = cat(1, RTs, secs - img_onset_t);

                    %% Case I: if 'h' is pressed
                    if keyCode(KbName('4'))
                        if any(j == [3, 8, 10])
                            old_responses = cat(1, old_responses, false);
                            correct_response = false;
                            misses = true;
                        else
                            old_responses = cat(1, old_responses, false);
                            correct_response = true;
                            correctRejections= true; 

                        end

                        %% Case II: if 'j' is pressed
                    elseif keyCode(KbName('5'))
                        if any(j == [3, 8, 10])
                            old_responses = cat(1, old_responses, true);
                            correct_response = true;
                            hits = true;
                        else
                            old_responses = cat(1, old_responses, true);
                            correct_response = false;
                            falseAlarms = true;
                        end
                    end
                end
                response_registered = true;
            else
                response_registered = false;
            end
        end

        img_offset_t = Screen('Flip', window, req_img_offsets_TR(j) + testrun_start_time - .5 * ifi);
        % Prepare drawing of dots already prior to experimental loop for ISI
        if j < n_practice_trials
            Screen('DrawDots', window, [xCenter; yCenter], 10, black, [], 2);
        end

        % Close image when it is no longer needed
        Screen('Close', texture);

        %image_onset_t relative to testrun starttime
        img_onset_t= img_onset_t - testrun_start_time;

        %image_onset_t relative to testrun starttime
        img_offset_t= img_offset_t - testrun_start_time;

        %Save image name to struct
        data.testRunImages(j).currentImageFilename = currentImageFilename;

        % Save image onset and offset times to struct
        data.testRunImages(j).img_onset = img_onset_t;
        data.testRunImages(j).img_offset = img_offset_t;

        % Save key responses, i.e., reaction times and recognition judgment
        data.testRunImages(j).RTs = RTs;
        data.testRunImages(j).old_responses = old_responses;
        % Saves whether person detected old and new image presetation
        % correctly - Hits and Correct recjection
        data.testRunImages(j).correct_response = correct_response;
        data.testRunImages(j).hits = hits;
        data.testRunImages(j).misses = misses;
        data.testRunImages(j).falseAlarms = falseAlarms;
        data.testRunImages(j).correctRejections = correctRejections;

        % Calculate the actual duration of stimulus presentation
        data.testRunImages(j).Image_actual_dur = img_offset_t - img_onset_t;

    end
    %Save data of correct and incorrect counts accordingly
    images_tab_TR = struct2table(data.testRunImages( ...
        1:n_practice_trials));

    % What was the percentage of correct responses in this block?


    %Increment the images tab hits
    %images_tab_TR.hits =  images_tab_TR.hits+1;

    % Set a reasonable feedback score to scale
    images_tab_TR.hits

    %testrun_crct_resp = mean([images_tab_TR.hits]/3)*5;
   testrun_crct_resp = mean([images_tab_TR.correct_response]) * 100;
   data.testrun_crct_resp = testrun_crct_resp; 

    %% Feedback text
    if data.currentLanguage == "English"

        % Provide feedback at the end
        feedbackMessage = sprintf('Your performance of detecting the old_images correctly in this block is rated with %.01f%%. Press any key to continue.', testrun_crct_resp);
        disp(feedbackMessage); % Display the feedback in the MATLAB command window
        % Display feedback to the participant on the screen
        DrawFormattedText(window, feedbackMessage, 'center', 'center', white);

    elseif data.currentLanguage == "German"
        % Provide feedback at the end with floating point number with one
        % decimal place
        feedbackMessage = sprintf('Sie haben in diesem Block %.01f%% der Fälle richtig beantwortet. Drücken Sie eine beliebige Taste, um fortzufahren.', testrun_crct_resp);
        disp(feedbackMessage); % Display the feedback in the MATLAB command window
        % Display feedback to the participant on the screen
        DrawFormattedText(window, feedbackMessage, 'center', 'center', white);
    end

    Screen('Flip', window);
    KbStrokeWait; % Wait for a key press to continue

    if data.currentLanguage == "English"
    % Display a message to proceed to the main experimental loop
    DrawFormattedText(window, 'Test Run Complete. Press any key to continue.', 'center', 'center', white);
    % Prepare drawing of dots already prior to experimental loop for ISI

    elseif data.currentLanguage == "German"
    % Display a message to proceed to the main experimental loop
    DrawFormattedText(window, 'Probedurchlauf beendet. Drücke eine beliebige Taste um fortzufahren. ', 'center', 'center', white);
    % Prepare drawing of dots already prior to experimental loop for ISI
    end

    Screen('Flip', window);
    KbStrokeWait;
else
    if data.currentLanguage == "English"
    % Display a message to proceed to the main experimental loop
    DrawFormattedText(window, 'Still feeling focused? Press any key to continue with the experiment.', 'center', 'center', white);
    % Prepare drawing of dots already prior to experimental loop for ISI

    elseif data.currentLanguage == "German"
    % Display a message to proceed to the main experimental loop
    DrawFormattedText(window, 'Noch fokussiert? Drücken Sie eine beliebige Taste, um fortzufahren.', 'center', 'center', white);
    % Prepare drawing of dots already prior to experimental loop for ISI
    end

    Screen('Flip', window);
    KbStrokeWait;
end

%-----------------------------------------------------------------------
%   Set number of trials (presented pictures 3s + 1s blank
%              space with fixation dot)
%-----------------------------------------------------------------------

% definition target size of resized image to fit given visual angle
% of 8.4
target_size = [data.pix_siz, data.pix_siz];

%----------------------------------------------------------------------
%      Create a Image Sequence Matrix (which images are shown in the
%      experiment with all blocks for one participant)
%----------------------------------------------------------------------


img_seq_path = sprintf('img_seq_sub%i_%iyrs_%s_%i.mat', data.nr, ...
    data.age, data.gender, todays_date);

% Do the retrieval of sequence matrix only in block 1
if data.currentBlockNumber==1

    % Call the  get_img_seq function from another file
    [img_seq, img_is_old] = get_img_seq(n_imgs, old_prob, n_blocks);
    %img_seq = get_img_seq(pic_ids, old_prob);


    % Complete image sequence (including image repetitions) requires this
    % many trials:
    n_trials = length(img_seq);

    %Debugging message that the division of n_trials (amount of trials) by
    %n_blocks(amount of blocks) should be equal to 0 --> if not true an error
    %message will be printed
    assert(mod(n_trials, n_blocks)==0, sprintf( ...
        ['You requested %i unique images with %.03f percent ' ...
        'probability for''old'' responses, which yielded %i\n' ...
        'trials in total. This is not integer-divisible by the ' ...
        'number of requested blocks, which is %i'], n_imgs, ...
        old_prob, n_trials, n_blocks));

    % change structre of img_seq into an structure arry with trials =n -->
    % n rows (each row has the field name pic) with value 1
    img_seq_matrix(n_trials).pic = [];

    for trial_num = 1:n_trials
        img_seq_matrix(trial_num).pic = img_seq(trial_num);
        img_seq_matrix(trial_num).img_path = pictures{img_seq(trial_num)};
        img_seq_matrix(trial_num).img_is_old = img_is_old(trial_num);
    end

    % Display the result or use it as needed
    disp('Generated Image Sequence:');
    disp(img_seq);

    save(img_seq_path, 'img_seq_matrix');
else
    load(img_seq_path, 'img_seq_matrix');
end

n_trials_per_block = length(img_seq_matrix) / n_blocks;

% Define adapted size of image with (8.4 degrees of visual angle like
% in Allen et al (2021)
trial_dur = requested_fixationDuration + requested_stimulusDuration;


req_img_onsets = (0:(n_trials_per_block-1)) * trial_dur + requested_fixationDuration;
req_img_offsets = req_img_onsets + requested_stimulusDuration;

% Text for starting calibration and validation

%%%%%%%%%%%%%%%%%%%%%%%% EYELINK Part III %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------------------------------------------
%            calibration and validation for each block
%------------------------------------------------------------------
if use_of_eyelink
    DrawFormattedText(window,'Press any key to continue with the calibration and validation', ...
        'center', 'center', [255 255 255]);
    fprintf('pre calibration: block %i \n', data.currentBlockNumber);

    EyelinkDoTrackerSetup(el);
    disp('calibration done')
end

if data.currentLanguage == "English"
    % Display a message to proceed to the main experimental loop
    DrawFormattedText(window, ['\n Ready for the next run: Remember please press "4" if you recognize a new image '...
    '\n and "5" if you recognize an old image. Press any key to continue.'], 'center', 'center', white);
    % Prepare drawing of dots already prior to experimental loop for ISI

elseif data.currentLanguage == "German"
    % Display a message to proceed to the main experimental loop
    DrawFormattedText(window, ['\n Sind Sie bereit für den nächsten Durchlauf? Drücken Sie bitte "4" wenn Sie ein neues Bild erkennen '...
    '\n und "5" wenn Sie ein altes Bild erkennen. Drücken Sie eine beliebige Taste um forzufahren.'], 'center', 'center', white);
    % Prepare drawing of dots already prior to experimental loop for ISI
end

Screen('Flip', window);
KbStrokeWait;

% if everything is set up
% PRESS "O" key (the letter) to start experiment -> output/ record
% on the host pc (eyetracker pc)

%%%%%%%%%%%%%%%%%%%%%%%% EYELINK Part IV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------------------------------------------
% Further eyelink code - Drift Correction and start of recording
%------------------------------------------------------------------
if use_of_eyelink
    % always put Host PC in offline mode
    % -> otherwise it might miss the StartRecording command!
    Eyelink('Command', 'set_idle_mode');
    WaitSecs(0.05);

    Eyelink('StartRecording');
    % record a few samples before we actually start displaying
    % otherwise you may lose a few msec of data
    WaitSecs(0.1);

    % get eye that's tracked
    % returns 0 (LEFT_EYE), 1 (RIGHT_EYE) or 2 (BINOCULAR)
    data.eye_used = Eyelink('EyeAvailable'); % FIXME
end
%% ------------------------------------------------------------------
%-----------------------------------------------------------------------
%                Start of experimental loop and eye tracking
%-----------------------------------------------------------------------

%%%%%%% Start of Experimental loop pictures an 3s On/Off trials %%%%%%%%%

% Query the frame duration
ifi = Screen('GetFlipInterval', window);

% prepare drawing dots in blank space (1s) = ISI
%%%%% preparation of drawing of dots and calculating target_size;

block_start_time = Screen('Flip', window);
Screen('DrawDots', window, [xCenter; yCenter], 10, black, [], 2);

% Save trial information in advance to data log
data.block_start_time = block_start_time;
data.req_img_onsets = req_img_onsets;
data.req_img_offsets = req_img_offsets;


block_first_trial_num = (data.currentBlockNumber - 1) * ...
    n_trials_per_block + 1;
block_last_trial_num = data.currentBlockNumber * n_trials_per_block;

trial_counter = 1;

for k = block_first_trial_num:block_last_trial_num

    %------------------------------------------------------------------
    %                     Get the timing
    %------------------------------------------------------------------

    data.images(k).isi_on = Screen('Flip', window, requested_fixationDuration)-block_start_time;

    % Get the presentation time of the ISI
    data.images(k).isi_off = WaitSecs(1)-block_start_time;
    % isRecDuration = complete time interval ISI was shown (supposed to
    % be 1s)
    data.images(k).isiRecDuration = data.images(k).isi_off - data.images(k).isi_on;
    %% Set up Calibration and Eyetracking
    % Now within the scope of each trial;

    %----------------------------------------------------------------------
    %            Determine the pictures presented in the trial
    %----------------------------------------------------------------------

    % Read picture
    %         curr_pic = char(img_seq_matrix(k).pic);
    cur_img_path = img_seq_matrix(k).img_path;

    % Current picture
    c_pic = imread(cur_img_path);

    %------------------------------------------------------------------
    %                    Define adapted size of image
    %------------------------------------------------------------------

    % Resize the image
    resized_image = imresize(c_pic, target_size);

    %------------------------------------------------------------------
    %                Show image on screen with make texture
    %------------------------------------------------------------------

    c_picTexture = Screen('MakeTexture', window, resized_image);

    Screen('DrawTexture', window, c_picTexture, [], [], 0);
    img_onset_t = Screen('Flip', window, data.req_img_onsets(trial_counter) + data.block_start_time - .5 * ifi);
    %% detect key presses during while loop and see if it's time to break the loop
    response_registered = false;
    old_responses = []; % Initialize to empty vector to be filled with key presses
    other_responses = []; % Initialize to empty vector to be filled with key presses

    % ResponseTime (RTs) caputures the moment the image was displayed ('image_onset_t) until a key press ('3' or '4' is detected)
    RTs = []; % Initialize to empty vector with response times relative to image onset

    % Setting correct response by default in false so that it covers also
    % those cases
    correct_response = false;
    hits = false;
    misses = false;
    falseAlarms = false; 
    correctRejections = false; 

    %%%%%%%%%%%%%%%%%%%%%%%% EXPERIMENTAL While Loop %%%%%%%%%%%%%%%%%%%%%%%%%
    % During Stimulus presentation (3s) if pressed '4' --> 'new image', if
    % pressed '5' --> 'old image' (should be pressed if subject thinks it has been presented to him before)
    % with pressing 'x' on the experimental keyboard the session can be aborte,
    % if any other key is pressed during stimulus presentation time the pressed
    % key(s) get stored in old_responses
    %% FIXME: prrof whether other key presses really get stored in old_responses

    while GetSecs() < data.block_start_time + data.req_img_offsets(trial_counter) - 1.5 * ifi

        [keyIsDown, secs, keyCode, ~] = KbCheck;
        if keyIsDown
            %%% abort experiment in pressing x during execution of while
            %%% loop
            if keyCode(KbName('x'))
                ShowCursor;
                sca;
                fclose all;
                Screen('CloseAll');
                return
            end

            %% FIXME: Implementation of adaption of feedback

            if (keyCode(KbName('4')) || keyCode(KbName('5'))) && response_registered==false
                RTs = cat(1, RTs, secs - img_onset_t);

                %% Case I: if img_is_old in img_seq_matrix is 0 --> picture is new

                if img_seq_matrix(k).img_is_old == 0

                    if keyCode(KbName('4'))
                        old_responses = cat(1, old_responses, false);
                        correct_response = true;
                        correctRejections = true;

                    elseif keyCode(KbName('5'))
                        old_responses = cat(1, old_responses, true);
                        correct_response = false;
                        falseAlarms = true;
                    end

                    %% Case II: if img_is_old in img_seq_matrix is 1 --> picture is old

                elseif img_seq_matrix(k).img_is_old == 1

                    if keyCode(KbName('4'))
                        old_responses = cat(1, old_responses, false);
                        correct_response = false;
                        misses = true;

                    elseif keyCode(KbName('5'))
                        old_responses = cat(1, old_responses, true);
                        correct_response = true;
                        hits = true; 
                    end
                end
            end

            response_registered = true;
        else
            response_registered = false;
        end

    end

    % Drawing of fixation dots again

    img_offset_t = Screen('Flip', window, (trial_counter) + block_start_time - .5 * ifi);
    if k < block_last_trial_num
        Screen('DrawDots', window, [xCenter; yCenter], 10, black, [], 2);
    end

    % Close image when it is no longer needed
    Screen('Close', c_picTexture);

    %------------------------------------------------------------------
    %               Data saving for experimental run
    %------------------------------------------------------------------

    % image_onset_t relative to the block starttime
    img_onset_t= img_onset_t - block_start_time;

    % image_onset_t relative to the block starttime
    img_offset_t= img_offset_t - block_start_time;

    % Save image to struct
    data.images(k).pictures = img_seq_matrix(k).img_path;

    % Save image onset and offset times to struct
    data.images(k).img_onset = img_onset_t;
    data.images(k).img_offset = img_offset_t;
    %________________________________________________________________________
    %% Save subject responses

    % Save key responses, i.e., reaction times and recognition judgment
    % data.images(k).response_registered = response_registered;

    % Save key responses, i.e., reaction times and recognition judgment
    data.images(k).RTs = RTs;

    % Saves whether person pressed '5', so detected an old response
    data.images(k).old_responses = old_responses;

    % Saves whether person detected old and new image presetation correctly
    data.images(k).correct_response = correct_response;
    data.images(k).hits = hits;
    data.images(k).misses = misses;
    data.images(k).correctRejections = correctRejections;
    data.images(k).falseAlarms = falseAlarms;

    % Saves whether person pressed any other key and if so, which one
    data.images(k).other_responses = other_responses;

    % Calculate the actual duration of stimulus presentation
    data.images(k).Image_actual_dur = data.images(k).img_offset  - data.images(k).img_onset;

    %% Increment trial counter
    trial_counter = trial_counter + 1;

end

%%%%%%%%%%%%%%%%%%%%%%%% EYELINK Part V %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if use_of_eyelink
    % adds 100 msec of data to catch final events
    WaitSecs(0.1);
    Eyelink('StopRecording');
end

%%%%%%%%%%%%%%%%%%%%%%%% EYELINK Part VII %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Eyetracking Code - Save Eye Tracking file (edf)
if use_of_eyelink
    Eyelink('CloseFile');

    % download the file
    try
        fprintf('Receiving data file ''%s''\n', edfFile );
        fileStatus = Eyelink('ReceiveFile');
        if fileStatus > 0
            fprintf('ReceiveFile status %d\n', fileStatus);
        end
        if 2==exist(edfFile, 'file')
            fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, pwd );
            % Diretly converting the file into an ascii format
        end

        % rdf is the exception variable
    catch rdf
        fprintf('Problem receiving data file ''%s''\n', edfFile );
        rdf;
    end
    % shut down Eyelink at end of Experiment
    % Eyelink('Shutdown');
    unixString = ['edf2asc ', edfFile ,'.edf'];
    system(unixString);
end


% Save data of correct and incorrect counts accordingly
images_tab = struct2table(data.images( ...
    block_first_trial_num:block_last_trial_num));

% What was the percentage of correct responses in this block?
block_crct_resp = mean(images_tab.correct_response) * 100;
data.block_cct_resp = block_crct_resp;

if data.currentLanguage == "English"

    % Provide feedback at the end
    feedbackMessage = sprintf('You answered correctly in %.01f%% of the trials in this block. Press any key to complete the block and start a break.', block_crct_resp);
    disp(feedbackMessage); % Display the feedback in the MATLAB command window

elseif data.currentLanguage == "German"

    % Provide feedback at the end
    feedbackMessage = sprintf('Sie haben in diesem Block %.01f%% der Fälle richtig beantwortet. Drücken Sie eine beliebige Taste, um den Block abzuschließen und eine Pause einzulegen.', block_crct_resp);
    disp(feedbackMessage); % Display the feedback in the MATLAB command window

end

% Display feedback to the participant on the screen
DrawFormattedText(window, feedbackMessage, 'center', 'center', white);
Screen('Flip', window);
KbStrokeWait; % Wait for a key press to continue

%%%%%%%%%%%%%%%% Show final ccreen of block and user feedback %%%%%%%%%%%%%%%%%%%%%%

%% TODO 8: Adapt to Blocknumbers n_blocks -1

% Draw all the text in one go

DrawFormattedText(window, text.endOfBlock ,'center', 'center', white);

% Flip to the screen
Screen('Flip', window);

% Wait for last key press to close the screen
KbStrokeWait;

% Close all open images
Screen('CloseAll');

% Close the screen
sca;

%_____________________________________________________________________________________________
%%%%%%%%%%%%%%%%%%%% Saving data to struct file (opens in matlab) %%%%%%%%%%%%%%%%%%%%%%%

%% According to trial data - What operating system do you use?
% Switch case determines what operating system is used --° Linux, macOS or
% Windows

% save data struct to a file with substructs testrun and images

switch computer ()
    %% Linux
    case 'GLNXA64'
        filename_data_mat= sprintf('/home/lab/Documents/lheinemann-ET-NSD-cocoImg/subjectData/data_subj_%i_%iyrs_%s_%i_%i.mat', data.nr, ...
            data.age, data.gender, todays_date, data.currentBlockNumber);
        save(filename_data_mat, "data");

        %% Apple (macOS)
    case 'MACI64'
        %         filename_data_mat= sprintf('/Volumes/bartels_data/lheinemann/subjectData/data_subj_%i_%iyrs_%s_%i_%i.mat', data.nr, ...
        %             data.age, data.gender, todays_date, data.currentBlockNumber);
        %         save(filename_data_mat, "data");

        %Code for saving the files locally on the macStudio in personal account

        filename_data_mat= sprintf('/Users/lisaheinemann/Documents/MATLAB/subjectData/data_subj_%i_%iyrs_%s_%i_%i.mat', data.nr, ...
            data.age, data.gender, todays_date, data.currentBlockNumber);
        save(filename_data_mat, "data");


        %% Windows
    case 'PCWIN64'
        filename_data_mat= sprintf('C:/Users/Lisa/MATLAB Drive/Matlab/subjectData/data_subj_%i_%iyrs_%s_%i_%i.mat', data.nr, ...
            data.age, data.gender, todays_date, data.currentBlockNumber);
        save(filename_data_mat, "data");
end

fprintf('\n Closing experiment screen and everything and Data is saved ...\n');

catch
    % If "try" cannot be executed, do this

%%%%%%%%%%%%%%%%%%%%%%% EYELINK Part VIII %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Close the eyetracker
    if use_of_eyelink
        Eyelink('Shutdown');
    end
%
    % Close all screens in case of an error
    Screen('CloseAll');
    % Close all open files
    fclose('all');
    psychrethrow(psychlasterror);
    fprintf('\nClosing error screen and everything due to an error...\n');
end

%% IMPORTANT: local function to receive img_seq
%------------------------------------------------------------------------
%  Local function for ramdomized calculated repetition with setted probability
%                          of old_responses
%------------------------------------------------------------------------

function [img_seq, img_old] = get_img_seq(n_imgs, old_prob, n_blocks)

%n_imgs = length(pic_ids);

%_______________________________________________________________________

%%%% 1st part

% How many images should be repeated?
n_imgs2rep = n_imgs / ((1 / old_prob) - 1);

assert(round(n_imgs2rep) == n_imgs2rep, 'Number of images and ''old'' probability are incompatible.');
assert(mod(n_imgs + n_imgs2rep, n_blocks) == 0, 'Cannot present the requested number of images in B blocks.');
% Maybe not neccessary? Proof
assert(mod(n_imgs, n_blocks) == 0, 'Number of images and number of blocks are incompatible.');

% Print information to debug
disp(n_imgs);
disp(n_imgs2rep);
disp(n_blocks);

% How many images will be presented per block?
n_imgs_per_block = (n_imgs + n_imgs2rep) / n_blocks;

% Print n_imgs_per_block to debug
disp(n_imgs_per_block);

% Which images should be repeated?
img2rep_ids = randperm(n_imgs, n_imgs2rep);

%_______________________________________________________________________

%%%% 2nd part

% The sequence consists of the following image IDs.
img_ids = cat(2, 1:n_imgs, img2rep_ids);

% First, distribute the repeated images uniformly
% across blocks: block_ids assigns each image to a
% single block. The assignment for the nonrepeated
% images is still unknown -- values are hence nan.
% The repeated images will be uniformly assigned
% to the B blocks.
block_ids = cat(2, nan(1, n_imgs), repmat(1:n_blocks, [1, n_imgs2rep / n_blocks]));

% The allocation of images looks as follows now:
% imagesc([img_ids / max(img_ids); block_ids]);

%_______________________________________________________________________

%%% 3rd part

candidate_assignments = kron(1:n_blocks, ones(1, n_imgs / n_blocks));

% Determine in what block the first presentation
% will occur
for seq_id = n_imgs+1:length(img_ids)
    % What is the ID of the image that corresponds
    % the current element of the sequence?
    cur_id = img_ids(seq_id);

    % In what block will the image be repeated?
    rep_block_id = block_ids(seq_id);

    % Need to present the first image in the same
    % block or earlier
    compatible_candidate_assignments = ...
        candidate_assignments(candidate_assignments ...
        <= rep_block_id);

    % Randomly choose a block where the first image
    % may appear
    cand_id = randi(length(compatible_candidate_assignments));
    block_ids(cur_id) = compatible_candidate_assignments(cand_id);

    % Update the vector of candidate assignments by
    % removing the selected candidate
    candidate_assignments = candidate_assignments( ...
        [1:cand_id - 1, cand_id + 1:end]);

end

% The assignment of images to blocks looks like this now:
% imagesc([img_ids; block_ids])

%_______________________________________________________________________

%%% 4th part

% The remaining images can be randomly assigned to
% fill up all blocks
for seq_id = find(isnan(block_ids))
    cur_id = img_ids(seq_id);

    % Randomly choose a block where the first image
    % may appear
    cand_id = randi(length(candidate_assignments));
    block_ids(cur_id) = candidate_assignments(cand_id);

    % Update the vector of candidate assignments by
    % removing the selected candidate
    candidate_assignments = candidate_assignments( ...
        [1:cand_id - 1, cand_id + 1:end]);

end

% The assignment of images to blocks looks like this now:
% imagesc([img_ids; block_ids]);

%_______________________________________________________________________

%%% 5th part

% Sort images depending on what block they appear in.
img_seq = [];
% n_blocks = 4
for blk = 1:n_blocks
    cur_block_seq = img_ids(block_ids==blk);
    cur_block_seq_shuffled = cur_block_seq(randperm(n_imgs_per_block));
    img_seq = [img_seq cur_block_seq_shuffled];
end


% Label each image based on whether it requires an "old"
% judgment or not.
img_old = zeros(1, n_imgs + n_imgs2rep);
for img_id = img_ids(n_imgs+1:end)
    img_ids_found = find(img_seq==img_id);
    assert(length(img_ids_found)==2, 'there must be two IDs');
    img_old(img_ids_found(2)) = true;
end

% Compute probability for "old" response for each block
% separately
old_prob_per_block = zeros(1, n_blocks);
for blk = 1:n_blocks
    block_start = (blk-1) * n_imgs_per_block + 1;
    block_end = blk * n_imgs_per_block;
    cur_old_prob = mean(img_old( ...
        block_start:block_end));
    old_prob_per_block(blk) = cur_old_prob;
    fprintf('Probability for ''old'' response in block %i is %.04f\n', blk, cur_old_prob);
end

%  Display the data in array C as an image that uses the full range of
%  colors in the colormap via imagesc([img_seq; img_old])
% The images will be presented in the following order
imagesc([img_seq; img_old])
