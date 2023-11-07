function Mu = friction_calc(Parms)
   a = Parms.a;
   n = Parms.n;
   dt = mean(diff(Parms.Time));
   k = Parms.k;
   
   x = linspace(-a, a, n+1);
   
   debug_time = 13.77;
   debug = -1;
   
   Prev = struct;
   Prev.idx = 0;
   Prev.q = zeros(size(x));
   Prev.qd = zeros(size(x));
   
   Mu = zeros(size(Parms.Time));
   
   wb = waitbar(0,'Calculating the coeff. of friction...');
   tic;
   for i = 1:length(Parms.Time)
       if (Parms.Time(i) >= debug_time) && (debug == 0)
           figure;
           fplot(@(x) minimize(x,Parms,Prev), [0, 2.5]);
           debug = 1;
       end
       
       %opts = optimoptions('fmincon', 'Display', 'off');
       mu_vect = fmincon(@(x) minimize(x,Parms,Prev), [2, 1.5], [], [], [],...
           [], [0, 0], [2.5, 2.5]);
       mu = mu_vect(1);
       mu_sl = mu_vect(2);
       
       q_next = Prev.q + Prev.qd*dt;
       qd_next = zeros(size(x));
       qcrp = 3*Parms.Fz(i)*mu/(4*a^3*k)*(a^2 - x.^2);
       q_slip = 3*Parms.Fz(i)*mu_sl/(4*a^3*k)*(a^2 - x.^2);
       omega = Parms.Omega(i);
       for ii = 1:n+1
            if abs(q_next(ii)) < qcrp(ii)
                qd_next(ii) = (Parms.l - x(ii))*omega;
            else
                qd_next(ii) = 0;
                q_next(ii) = (q_slip(ii) - Parms.epsilon)*sign(q_next(ii));
            end
       end
       
       Mu(i) = mu;
       Prev.idx = i;
       Prev.q = q_next;
       Prev.qd = qd_next;
       waitbar(i/length(Parms.Time),wb)
   end
   calcTime = toc
   close(wb);
   
   %% Internal function
    function cost = minimize(x, Parms, Prev)
        mu = x(1);
        mu_slip = x(2);
        
        if mu_slip > mu
            cost = 100;
            return;
        end
        
        idx = Prev.idx + 1;
        q_prev = Prev.q;
        qd_prev = Prev.qd;
        fz = Parms.Fz(idx);
        a = Parms.a;
        k = Parms.k;
        n = Parms.n;
        x = linspace(-a, a, n+1);
        omega = Parms.Omega(idx);
        
        qcrp = 3*fz*mu/(4*a^3*k)*(a^2 - x.^2);
        q_slip = 3*fz*mu_slip/(4*a^3*k)*(a^2 - x.^2);
        q = q_prev + qd_prev*dt; % Euler formula --> simple, but quick
        qd = zeros(1,n+1);
        
        for j = 1:n+1
            if abs(q(j)) < qcrp(j)
                qd(j) = (Parms.l - x(j))*omega;
            else
                qd(j) = 0;
                q(j) = (q_slip(j) - Parms.epsilon)*sign(q(j));
            end
        end
        
        Mz = 0;
        for j = 2:n+1
            Mz = Mz + (Parms.l - (x(j) + x(j-1))/2)*(q(j) + q(j-1))/2*Parms.dx;
        end
        Mz = -k*Mz;
        
        cost = (Mz - Parms.Mz(idx))^2;
    end
end