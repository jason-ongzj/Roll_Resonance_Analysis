% Class written for RollResonance as a way to encapsulate data. This 
% is so that functions do not need to be written with too many input
% arguments required for computation.

classdef RollResonance < handle
    properties
        % Astos input data
        Time, Ixx, Iyy, Pss, V, mach, q, T, m, CoM, CoP, Cd, rho
        AstosSS,AstosAoA
        
        % Reference quantities
        D, Sref, cant
        
        % Assumed quantities
        eps, gamma, Toffset, CGoffset, staticAoA, mu, lambda
        
        % Interpolated quantities
        CL_Delta, CL_P, CN_A, Mach_CLP_CLDelta, Mach_CNA
        CMQ_x, CMQ_y, CMQ_z
        
        % Variables for calculation of roll and pitch frequencies
        P, Pdot, roll, phase, trimAoA, POmegaRatio, omega, Prk4, PssCalc
        
        % Variables for roots of linear 2nd ODE and trimArm
        rootAReal, rootAImag, rootBReal, rootBImag, trimArmReal, trimArmImag
        
        % Variables for calculation of K1, K2 and K3 arms
        sideslipRate, sideslip, AoARate, AoA
        K1real, K1imag, K2real, K2imag, K3real, K3imag
    end
    
    methods
        % Initialize class object with interpolation of aerodynamic
        % characteristics.
        % Initialize object with: RollResonance("CN_A.txt", "CL_Delta.txt", "CL_P.txt")
        function self = RollResonance(cna_file, cl_delta_file, clp_file, cmq_file)
            cl_delta_data = importdata(cl_delta_file, "\t", 1);
            cl_p_data = importdata(clp_file, "\t", 1);
            cna_data = importdata(cna_file, "\t", 1);
            cmq_data = importdata(cmq_file, "\t", 2);
            
            self.CL_Delta = cl_delta_data.data(:, 2);
            self.CL_P = cl_p_data.data(:, 2);
            self.Mach_CLP_CLDelta = cl_p_data.data(:, 1);
            
            self.CN_A = cna_data.data(:, 2)*-1;
            self.Mach_CNA = cna_data.data(:, 1);
            
            % Interpolation for CMQ
            % Interpolate using interp2(self.CMQ_x, self.CMQ_y, self.CMQ_z, CG_from_aft, mach)
            % E.g. interp2(self.CMQ_x, self.CMQ_y, self.CMQ_z, 1.7555, 0.5)
            
            fid=fopen(cmq_file);
            % CG from aft
            self.CMQ_x = str2double(split(fgetl(fid)));
            % Mach number
            self.CMQ_y = str2double(split(fgetl(fid)));
            fclose(fid);
            self.CMQ_z = cmq_data.data;
        end
        
        function InitInputData(obj, astos_output)
            file_input = importdata(astos_output, '\t', 1);
            obj.Time = file_input.data(:, 1);
            obj.Ixx = file_input.data(:, 2);
            obj.Iyy = file_input.data(:, 3);
            obj.Pss = file_input.data(:, 4);
            obj.V = file_input.data(:, 6)*1000;
            obj.mach = file_input.data(:, 7);
            obj.q = file_input.data(:, 8);
            obj.T = file_input.data(:, 9)*1000;
            obj.m = file_input.data(:, 10)*1000;
            obj.CoM = file_input.data(:, 11);
            obj.CoP = file_input.data(:, 12);
            obj.Cd = file_input.data(:, 13);
            obj.rho = file_input.data(:, 18);
            
            % Solution for epicyclic theory within Astos
            obj.AstosSS = file_input.data(:, 14);
            obj.AstosAoA = file_input.data(:, 16);
        end
        
        function InitializeVars(obj)
            % Roll and pitch frequencies
            obj.P = zeros(length(obj.Time), 1);
            obj.Prk4 = zeros(length(obj.Time), 1);
            obj.PssCalc = zeros(length(obj.Time), 1);
            obj.Pdot = zeros(length(obj.Time), 1);
            obj.roll = zeros(length(obj.Time), 1);
            obj.phase = zeros(length(obj.Time), 1);
            obj.trimAoA = zeros(length(obj.Time), 1);
            obj.POmegaRatio = zeros(length(obj.Time), 1);
            obj.omega = zeros(length(obj.Time), 1);
            
            % Roots of linear 2nd ODE and trimArm
            obj.rootAReal = zeros(length(obj.Time), 1);
            obj.rootAImag = zeros(length(obj.Time), 1);
            obj.rootBReal = zeros(length(obj.Time), 1);
            obj.rootBImag = zeros(length(obj.Time), 1);
            obj.trimArmReal = zeros(length(obj.Time), 1);
            obj.trimArmImag = zeros(length(obj.Time), 1);
            
            % K1, K2 and K3 arms
            obj.sideslipRate = zeros(length(obj.Time), 1);
            obj.sideslip = zeros(length(obj.Time), 1);
            obj.AoARate = zeros(length(obj.Time), 1);
            obj.AoA = zeros(length(obj.Time), 1);
            
            obj.K1real = zeros(length(obj.Time), 1);
            obj.K1imag = zeros(length(obj.Time), 1);
            obj.K2real = zeros(length(obj.Time), 1);
            obj.K2imag = zeros(length(obj.Time), 1);
            obj.K3real = zeros(length(obj.Time), 1);
            obj.K3imag = zeros(length(obj.Time), 1);
        end
        
        function SetRefQuantities(obj, D, Sref, cant)
            obj.D = D;
            obj.Sref = Sref;
            obj.cant = cant;
        end
        
        function AssumeQuantities(obj, eps, gamma, Toffset, CGoffset, staticAoA, mu, lambda)
            obj.eps = eps;
            obj.gamma = gamma;
            obj.Toffset = Toffset;
            obj.CGoffset = CGoffset;
            obj.staticAoA = staticAoA;
            obj.mu = mu;
            obj.lambda = lambda;
        end
        
        function FindNaturalFrequency(obj)
            for i = 1:length(obj.Time)
                cna = interp1(obj.Mach_CNA, obj.CN_A, obj.mach(i), 'makima');
                obj.omega(i) = sqrt(-cna * (obj.CoM(i) - obj.CoP(i))...
                    * obj.q(i) * obj.Sref / obj.Iyy(i));
            end
        end
        
        function cmq = CMQInterpolation(obj, i)
           mach = obj.mach(i);
           if (mach < 0.5)
               cmq = interp2(obj.CMQ_x, obj.CMQ_y, obj.CMQ_z, obj.CoM(i), 0.5);
           elseif (mach > 4)
               cmq = interp2(obj.CMQ_x, obj.CMQ_y, obj.CMQ_z, obj.CoM(i), 4);
           else
               cmq = interp2(obj.CMQ_x, obj.CMQ_y, obj.CMQ_z, obj.CoM(i), mach, 'linear');
           end
           return;
        end
        
        function cna = CNAInterp(obj, i)
           mach = obj.mach(i);
           if (mach < 0.5)
               cna = interp1(obj.Mach_CNA, obj.CN_A, 0.5, 'linear');
           else
               cna = interp1(obj.Mach_CNA, obj.CN_A, mach, 'linear');
           end
           return;
        end
        
        function clp = CLPInterp(obj, i)
           mach = obj.mach(i);
           if (mach < 0.5)
               clp = interp1(obj.Mach_CLP_CLDelta, obj.CL_P, 0.5, 'linear');
           else
               clp = interp1(obj.Mach_CLP_CLDelta, obj.CL_P, mach, 'linear');
           end
           return;
        end
        
        function cldelta = CLDeltaInterp(obj, i)
           mach = obj.mach(i);
           if (mach < 0.5)
               cldelta = interp1(obj.Mach_CLP_CLDelta, obj.CL_Delta, 0.5, 'linear');
           else
               cldelta = interp1(obj.Mach_CLP_CLDelta, obj.CL_Delta, mach, 'linear');
           end
           return;
        end
        
        function b = CalcBVal(obj, i)
           m = obj.m(i);
           Ixx = obj.Ixx(i);
           Iyy = obj.Iyy(i);
           d = obj.D;
           mach = obj.mach(i);
           Sref = obj.Sref;
           q = obj.q(i);
           cna = interp1(obj.Mach_CNA, obj.CN_A, mach, 'linear');
           b = q*Sref*(-cna*(1 - Ixx/Iyy) - (m*d^2) * obj.CMQInterpolation(i)/Iyy) / ...
               (m * obj.V(i) * obj.omega(i));
           return;
        end
        
        function CalcTrimAoA(obj, i)
            P = obj.P(i)/360;
            w = obj.omega(i);
            b = obj.CalcBVal(i);
            denom = sqrt((1 - (P^2/w^2)*(1 - obj.Ixx(i)/obj.Iyy(i)))^2 + (b*P/w)^2);
            obj.trimAoA(i) = obj.staticAoA/denom;
        end
        
        function CalcPhaseAngle(obj, i)
           P = obj.P(i)/360;
           w = obj.omega(i);
           b = obj.CalcBVal(i);
           denom = 1  - (P^2/w^2)*(1-obj.Ixx(i)/obj.Iyy(i));
           obj.phase(i) = rad2deg(atan((b*P/w)/denom));
           obj.POmegaRatio(i) = P/w;
        end
        
        function CalcDynamicRoll(obj, i, start)
           d = obj.D;
           q = obj.q(i);
           Sref = obj.Sref;
           mach = obj.mach(i);
           Ixx = obj.Ixx(i);
           cant = obj.cant;
           cg_offset = obj.CGoffset;
           delta_t = obj.Time(i+1) - obj.Time(i);
           v = obj.V(i);
           trimAoA = obj.trimAoA(i);
           roll = obj.roll(i);
           gamma = obj.gamma;
           
           obj.CalcTrimAoA(i);
           obj.CalcPhaseAngle(i);
           
           cna = obj.CNAInterp(i);
           clp = obj.CLPInterp(i);
           cldelta = obj.CLDeltaInterp(i);
            
           A = q*Sref*d*clp*d*0.5/(v*Ixx);
           B = q*Sref*d*(cldelta*cant- cg_offset*cna*trimAoA*sin(deg2rad(roll - gamma)))/Ixx;
           
           if (obj.Time(i+1) <= 20.1)
               obj.P(i+1) = -B/A;
           else
               obj.P(i+1) = -B/A + (obj.P(2100) + B/A)*exp(A*(obj.Time(i)-obj.Time(2100)));
           end
           obj.roll(i+1) = obj.roll(i) + 0.5*(obj.P(i+1) + obj.P(i))*delta_t;
        end

        function CalcSteadyStateRoll(obj, i)
            cant = obj.cant;
            V = obj.V(i);
            d = obj.D;
            delta_t = obj.Time(i+1) - obj.Time(i);
            clp = obj.CLPInterp(i);
            cldelta = obj.CLDeltaInterp(i);
            obj.PssCalc(i) = -cldelta*cant*2*V/(clp*d);
            obj.roll(i+1) = obj.roll(i) + 0.5*(obj.P(i+1) + obj.P(i))*delta_t;
        end
        
%       Find sideslip and AoA values to determine coning motion of rocket
        function A1 = CalcA1(obj, i)
            m = obj.m(i);
            d = obj.D;
            Sref = obj.Sref;
            V = obj.V(i);
            I = obj.Iyy(i);
            g = 9.81;
            Isp = 233;
            T = obj.T(i);
            Toffset = obj.Toffset;
            q = obj.q(i);
            cna = obj.CNAInterp(i);
            cmq = obj.CMQInterpolation(i);
            A1 = q*Sref*(cna - m*d^2 * cmq/I)/(m*V) + T/(m*V) + ...
                    T*Toffset^2/(I*g*Isp);
            return;
        end
        
        function A2 = CalcA2(obj, i)
           P = obj.P(i)/360;
           Ixx = obj.Ixx(i);
           I = obj.Iyy(i);
           A2 = P * (2 - Ixx/I);
           return;
        end
        
        function B1 = CalcB1(obj, i)
            P = obj.P(i)/360;
            Sref = obj.Sref;
            d = obj.D;
            Ixx = obj.Ixx(i);
            I = obj.Iyy(i);
            cna = obj.CNAInterp(i);
            q = obj.q(i);
            x_margin = obj.CoM(i) - obj.CoP(i);
            B1 = cna*x_margin*q*Sref/I + P^2 * (1 - Ixx/I);
            return;
        end
        
        function B2 = CalcB2(obj, i)
            P = obj.P(i)/360;
            Pdot = obj.Pdot(i)/360;
            Ixx = obj.Ixx(i);
            I = obj.Iyy(i);
            m = obj.m(i);
            V = obj.V(i);
            Sref = obj.Sref;
            d = obj.D;
            g = 9.81;
            Isp = 233;
            T = obj.T(i);
            q = obj.q(i);
            Toffset = obj.Toffset;
            cna = obj.CNAInterp(i);
            cmq = obj.CMQInterpolation(i);
            B2 = -P*(q*Sref/(m*V))*(cna*(1-Ixx/I) - m*d^2/I*cmq) - ...
                P*(T/(m*V))*(1-Ixx/I) - P*T*Toffset^2/(I*g*Isp) - Pdot;
            return;
        end
        
        function C1 = CalcC1(obj, i)
            Sref = obj.Sref;
            gamma = deg2rad(obj.gamma);
            mu = deg2rad(obj.mu);
            d = obj.D;
            I = obj.Iyy(i);
            CGoffset = obj.CGoffset;
            Toffset = obj.Toffset;
            T = obj.T(i);
            eps = obj.eps;
            q = obj.q(i);
            cd = obj.Cd(i);
            cna = obj.CNAInterp(i);
            cmo = -cna * (obj.CoM(i) - obj.CoP(i)) * obj.staticAoA/d;
            C1 = q*Sref*d*(cmo*cos(obj.lambda) - cd*CGoffset*cos(gamma))/I + ...
                T*(CGoffset + Toffset)*eps*cos(mu)/I;
            return;
        end
        
        function C2 = CalcC2(obj, i)
            Sref = obj.Sref;
            gamma = deg2rad(obj.gamma);
            mu = deg2rad(obj.mu);
            d = obj.D;
            I = obj.Iyy(i);
            cna = obj.CNAInterp(i);
            cmo = -cna * (obj.CoM(i) - obj.CoP(i)) * obj.staticAoA/d;
            CGoffset = obj.CGoffset;
            Toffset = obj.Toffset;
            eps = obj.eps;
            T = obj.T(i);
            q = obj.q(i);
            cd = obj.Cd(i);
            C2 = q*Sref*d*(cmo*sin(obj.lambda) - cd*CGoffset*sin(gamma))/I + ...
                T*(CGoffset + Toffset)*eps*sin(mu)/I;
            return;
        end
        
        function CalcRoots(obj, i)
           A1 = obj.CalcA1(i);
           A2 = obj.CalcA2(i);
           B1 = obj.CalcB1(i);
           B2 = obj.CalcB2(i);
           C1 = obj.CalcC1(i);
           C2 = obj.CalcC2(i);
           
           A = A1 + A2*1j;
           B = -(B1 + B2*1j);
           C = C1 + C2*1j;
           
           % Roots
           a = -0.5*A + sqrt(A^2/4 - B);
           b = -0.5*A - sqrt(A^2/4 - B);
           
           obj.rootAReal(i) = real(a);
           obj.rootAImag(i) = imag(a);
           obj.rootBReal(i) = real(b);
           obj.rootBImag(i) = imag(b);
           
           trimArm = -C/B;
           obj.trimArmReal(i) = real(trimArm);
           obj.trimArmImag(i) = imag(trimArm);
        end
        
        function InitializeBC(obj, i, sideslip_rate, sideslip, aoa_rate, aoa)
           obj.sideslipRate(i) = sideslip_rate;
           obj.sideslip(i) = sideslip;
           obj.AoARate(i) = aoa_rate;
           obj.AoA(i) = aoa;
        end
        
        % Refer to Harold Vaughn's paper for tricyclic theory breakdown
        function CalcComplexPlane(obj, i)
           xi_rate_prev = obj.sideslipRate(i-1) + obj.AoARate(i-1)*1j;
           xi_prev = obj.sideslip(i-1) + obj.AoA(i-1)*1j;
           
           rootA = obj.rootAReal(i) + obj.rootAImag(i)*1j;
           rootB = obj.rootBReal(i) + obj.rootBImag(i)*1j;
           trimArm = obj.trimArmReal(i) + obj.trimArmImag(i)*1j;
           
           p = obj.P(i)/360;
           delta_t = obj.Time(i) - obj.Time(i-1);
           
           % Compute K1, K2  and K3 arms
           K3 = trimArm;%/exp(1j*p*delta_t);
           K1 = (xi_rate_prev - rootB*xi_prev - (1j*p - rootB)*K3)/(rootA - rootB);
           K2 = (xi_rate_prev - rootA*xi_prev - (1j*p - rootA)*K3)/(rootB - rootA);
           
           obj.K1real(i) = real(K1);
           obj.K1imag(i) = imag(K1);
           obj.K2real(i) = real(K2);
           obj.K2imag(i) = imag(K2);
           obj.K3real(i) = real(K3);
           obj.K3imag(i) = imag(K3);
           
           xi_current = K1*exp(rootA*delta_t) + K2*exp(rootB*delta_t) + K3*exp(1i*p*delta_t);

           obj.AoA(i) = imag(xi_current);
           obj.sideslip(i) = real(xi_current);
           
           % Compute sideslip and AoA rates
           xi_dot = K1*rootA*exp(rootA*delta_t) + K2*rootB*exp(rootB*delta_t) + K3*1i*p*exp(1i*p*delta_t);
           obj.AoARate(i) = imag(xi_dot);
           obj.sideslipRate(i) = real(xi_dot);
           
           % Update Pdot
           obj.Pdot(i) = (obj.P(i)-obj.P(i-1))/delta_t;
        end
    end
end
    