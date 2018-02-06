%MEAN PHASE COHERENCE
%CODED BY: PUNEET DHEER(pd)
%Date: 05-08-2017
%
%References:
%[1]Mean phase coherence as a measure for phase synchronization and its
%   application to the EEG of epilepsy patients. (by Florian Mormann,2000)
%
%INPUT:
% signal = column wise signal
% Ws = Window size in sample point
% Lp = Leave point or shifted by
% demean = 'yes' or 'no'
%
%OUTPUT:
% w_m_p_c = window wise mean phase coherence in matrix form
% channel_pairs = channel pair wise calculation

%%
function [w_m_p_c,channel_pairs]=MeanPhaseCoherence(signal,Ws,Lp,demean)
tic

Lw=1;
Z=Ws;
X=signal; 
no_of_channels=size(X,2);

windows=ceil((length(X)-Ws+1)/Lp); 

    for aa=1:windows

        windowed_data=X(Lw:Z,:);
        Lw=Lw+Lp;
        Z=Z+Lp;
        % Mean subtraction (demean), use only when it is required
        if strcmp('yes',demean)
            windowed_data = windowed_data-(ones(length(windowed_data), 1)*mean(windowed_data));
        end
        
        Phase_ts = hilbert_ts(windowed_data);
        rawMPC = abs(MPC(Phase_ts));
%         m_p_c{1,aa} = rawMPC-diag(diag(rawMPC)); %to make diagonal zero (optional)
        w_m_p_c{1,aa} = rawMPC;
        
        figure(1)
        title('MEAN PHASE COHERENCE')
        pcolor(w_m_p_c{1,aa})
        %shading interp; 
        shading flat; 
        axis ij
        colormap(jet)
                   
%         pause(0.1)
        

    end
    
   channel_pairs=[]; 
    for i=1:size(w_m_p_c,2)
        w_data=w_m_p_c{1,i};
        ch=1;
        for j=1:no_of_channels
            for k=j+1:no_of_channels
                c_p(ch,:)=w_data(j,k);
                ch=ch+1;
            end
        end
        channel_pairs=[channel_pairs c_p];
    end
    
    
 toc
end


function Ph_ts = hilbert_ts(d_data)
    hts = hilbert(d_data); %hibert transformation (analytical signal)
    Ph_ts = atan2(imag(hts), real(hts)); %must be real
%     Ph_ts = angle(hts);
end


function windowed_mpc = MPC(Ph)
    [N, D] = size(Ph);
    %Ph = mod(Ph,2*pi); % 0 to 2pi
    windowed_mpc = ones(D,D);
    for dim1= 1:D
        for dim2= dim1+1:D
            rel_phase_diff=(Ph(:, dim1)-Ph(:, dim2)); %-pi to +pi
            windowed_mpc(dim1, dim2)= abs(mean(exp(1i* mod(rel_phase_diff,2*pi)))); % converted to[0,2pi)
%           windowed_mpc(dim2, dim1)= conj(windowed_mpc(dim1, dim2));
            windowed_mpc(dim2, dim1)= windowed_mpc(dim1, dim2);
        end
    end
end

