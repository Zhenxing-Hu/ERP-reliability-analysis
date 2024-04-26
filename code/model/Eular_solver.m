function result = Eular_solver(data,A,h,C,noise,noise_var,jitter,u)
% A:         system transition matrix
% h:         Iteration step
% C:         input strength
% noise:     pink noise
% noise_var: noise strength
% jitter:    subject-level latency jitter
% u:         system input
        for t=1:1999
            data(t+1,:) = data(t,:)+h*(data(t,:)*A+noise(t)*noise_var);
            if t >jitter+500
                data(t+1,:) = data(t+1,:)+ C*[u(t-jitter-500),0];%
            end
        end
        result = data;
end