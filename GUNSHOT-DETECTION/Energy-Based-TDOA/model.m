% create vector frames (512x1) and put in workspace
fs=48000;
frame_len = 4096;
t = (0:frame_len-1)'/fs;
pulse = zeros(frame_len,1); pulse(50)=1; pulse = conv(pulse,hann(20),'same'); pulse = pulse / max(abs(pulse));
audio1_vec = 0.2 * pulse;
audio2_vec = 0.5 * circshift(pulse,7);
audio3_vec = 0.7 * circshift(pulse,11);

assignin('base','audio1_vec',audio1_vec);
assignin('base','audio2_vec',audio2_vec);
assignin('base','audio3_vec',audio3_vec);
mic_positions = [0 2 1; 0 0 1.5];
c = 343;
energy_th = 0.05;
assignin('base','energy_th',energy_th);
fs = 48000; frame_len = 4096;
assignin('base','fs',fs); assignin('base','frame_len',frame_len);
set_param(bdroot,'Solver','FixedStepDiscrete','FixedStep',num2str(1/fs),'StopTime','1.0');
snippet_len = 512;
assignin('base','snippet_len',double(snippet_len));
