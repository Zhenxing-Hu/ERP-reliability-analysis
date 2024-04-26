
function [ERP_ses1,ERP_ses2,grand_average] = Simulation(param,varargin)
%%%%%%%%%%%%%%%%%%%%%%% Model description %%%%%%%%%%%%%%%%%%%%%%% 
%  Goal statement
%  In short, the objective of this simulation is to explain what we observed in real ERP data 
%  (i.e., the inconsistency between group effects and individual reliability) by biophysical neuron models. 
%  Such differential equation modeling technique has been commonly applied to studying the mechanism of EEG\ERP
%  generation (Hebbink et al., 2020)(Jansen and Rit, 1995)(Huang et al., 2011)(David et al., 2005)(Romagnoni et al., 2020). 
%  The design of this model does not aim to simulate the behavior of brain wave oscillation as realistically as possible, 
%  but to explain our observations with the simplest possible model.

%  Model assumptions:
%  Two main assumptions with justifications are listed below:

%  Assumptions #1: Only subject-level transmission delay for stimuli is introduced, and inter-trial 
%  variability is modeled as perceptual variability of input strength.

%  Assumptions #2: Based on the observation of the real data, there are no systematic state variances 
%  (i.e., between-session variances). Hence proportion of Var(State) is negligible.

%  Model specification: x' = Ax(t)+C*u(t)+e(t)
%  A is the state-transition matrix,The elements of A are selected empirically to satisfy 
%  the following loose requirements: 
%  1.	The system would converge to baseline after a while of the stimulus.
%  2.   The time scales of the simulated damped oscillations should be comparable with ERP waveform of real data.

%  The inputs can be divided into two parts: one is e(t), the EEG background activity simulated by 1/f noise 
%  (i.e., pink noise), u(t) refers to sensory stimulation as shown from equation (7) in our submitted manuscripts. 

%  Limitation:
%  Assumptions #3: For simplicity, system transition matrix A is fixed across different subjects and trials, which 
%  is questionable because realistically individual-specific parameterization should be adopted, as stated by (Seghier and Price, 2018).

% References:

%   David, O., Harrison, L., Friston, K.J., 2005. Modelling event-related responses in the brain. Neuroimage 25, 756每770. https://doi.org/10.1016/j.neuroimage.2004.12.030
%   Hebbink, J., van Gils, S.A., Meijer, H.G.E., 2020. On analysis of inputs triggering large nonlinear neural responses Slow-fast dynamics in the Wendling neural mass model: Slow-fast dynamics in the Wendling neural mass model. Commun. Nonlinear Sci. Numer. Simul. 83, 105103. https://doi.org/10.1016/j.cnsns.2019.105103
%   Huang, G., Zhang, D., Meng, J., Zhu, X., 2011. Interactions between two neural populations: A mechanism of chaos and oscillation in neural mass model. Neurocomputing 74, 1026每1034. https://doi.org/10.1016/j.neucom.2010.11.019
%   Jansen, B.H., Rit, V.G., 1995. Electroencephalogram and visual evoked potential generation in a mathematical model of coupled cortical columns. Biol. Cybern. 73, 357每366. https://doi.org/10.1007/BF00199471
%   Romagnoni, A., Colonnese, M.T., Touboul, J.D., Gutkin, B.S., 2020. Progressive alignment of inhibitory and excitatory delay may drive a rapid developmental switch in cortical network dynamics. J. Neurophysiol. 123, 1583每1599. https://doi.org/10.1152/jn.00402.2019
%   Seghier, M.L., Price, C.J., 2018. Interpreting and Utilising Intersubject Variability in Brain Function. Trends Cogn. Sci. 22, 517每530. https://doi.org/10.1016/j.tics.2018.03.003

    if ~isempty(varargin)
        param_range = param.(varargin{:});
    else
        for_range = 1; 
    end
    
    grand_average = zeros(for_range,1500);
    
    % Parameter initialization
    C_sub_mean = param.C_sub_mean;
    C_sub_var  = param.C_sub_var;
    C_trial_mean = param.C_trial_mean;
    C_trial_var = param.C_trial_var;
    jitter_sub_mean = param.jitter_sub_mean;
    jitter_sub_var = param.jitter_sub_var;
    noise_var = param.noise_var;
    nsub = param.nsub;
    ntrial = param.ntrial;
    h = param.iter_step;
     % State-transition matrix
    c = param.c; d = param.d;
    A = [c,d;-d,c];
    
    % System input
    t=(1:2000)/1000;
    a = 20;b=50;
    u = a.*t.*exp(-b*t);
    
     % Pink noise
    weight=ones(1999,1);
    f=linspace(0,1000,1999);
    f(1999:-1:1001)=f(2:1000);
    weight(2:1999)=1./f(2:1999);
    % Store simulated result
    result = zeros(for_range,1500,nsub,ntrial);
    
    for idx_range = 1:for_range 
        if for_range > 1 
            eval([varargin{:},' = param_range(idx_range)']);
            if contains(['c','d'],varargin{:})
                A = [c,d;-d,c];
            end
        end
        for sub=1:nsub
            clc;disp(['nsub:',num2str(sub),'   idx_range:',num2str(idx_range)]);
            C_sub = C_sub_mean+randn(1,1)*C_sub_var;
            jitter = jitter_sub_mean+floor((rand(1,1)-0.5)*jitter_sub_var);
            parfor trial=1:ntrial
                result_init = zeros(2000,2);
                C_trial = C_trial_mean + C_trial_var*randn(1,1);
                C = C_sub + C_trial;
                noise = real(ifft(fft(randn(1999,1)).*weight));
                result_temp =  Eular_solver(result_init,A,h,C,noise,noise_var,jitter,u);
                % baseline correction
                result_temp(:,1)=result_temp(:,1)-mean(result_temp(500+(1:500),1));
                result(idx_range,:,sub,trial)=result_temp(501:2000,1);
            end
        end
        ERP_ses1 =  squeeze(mean(result(:,:,:,1:floor(ntrial/2)),4));
        ERP_ses2 =  squeeze(mean(result(:,:,:,(1:floor(ntrial/2))+floor(ntrial/2)),4));
        grand_average(idx_range,:,:) = mean(squeeze(result(idx_range,:,:)),2);
    end
end


