function mvt_pdf = movementModel(method, varargin)

% compit 
speed_x = (0:1000)';

switch method
    case "gam"
        if nargin>2
            speed_y = gampdf(speed_x,varargin{1},varargin{2});
        else
            speed_y = gampdf(speed_x,7,7);
        end
        
    case "norm"
        speed_y = normpdf(speed_x,30, 10);

    case "energy"
        if nargin<3
            %         name                              mass       wing
            %         ____________________________    _______    _______ 
            %         'Willow Warbler'                 0.009      0.195  
            %         'Tree Pipit'                    0.0225       0.26  
            %         'Common Chiffchaff'             0.0075       0.18  
            %         'Spotted Flycatcher'             0.016       0.24  
            %         'Garden Warbler'                0.0195       0.22  
            %         'Common Whitethroat'             0.015       0.21  
            %         'European Pied Flycatcher'       0.012       0.23  
            %         'Common Redstart'                0.016      0.225  
            %         'Wood Warbler'                  0.0095       0.22  
            %         'Eurasian Blackcap'              0.017       0.23  
            %         'Song Thrush'                    0.083       0.34 
                
             m = .083; %[kg] mass of bird.
             B = .34; %[m] Wing span.
        else
            m = varargin{1}/1000; %[kg] mass of bird.
            B = varargin{2}/100; %[m] Wing span.
        end

        % Parameter for energy speding
        k=1.2; % [-] Induced power factor (p. 45).
        gcst = 9.81; % [ms-2] gravity constant
        airdens = 1; % Air density 
        Sb = @(m) 0.00813*m^0.666; % [m2] body frontal area
        CDb = 0.1; % [-] body drag coefficient (p. 51).
        Cpro = 8.4;
        Ra = 7;

        % Vt [m/s] true speed (bird speed-wind speed)
        f_eng = @(Vt) (2*k*(m*gcst)^2)./(Vt*pi*B.^2.*airdens) + airdens.*Vt.^3*Sb(m)*CDb/2 + Cpro/Ra*1.05*k^(3/4)*m^(3/2)*gcst^(3/2)*Sb(m)^(1/4)*CDb^(1/4)./airdens.^(1/2)/B^(3/2);

        speed_y = (1./f_eng(speed_x)).^.5;
end

% 
speed_y(speed_y<.005&speed_x<40)=.005;


%% Create pdf
% slow version
% mvt_pdf = @(x) interp1(speed_x,speed_y,x,'linear',eps);

% fasterer version with a rounding of speed to 1 km/h
mvt_pdf = @(x) speed_y(min(round(x)+1,numel(speed_y)));


%% Figure
% figure; plot(speed_x*1000/60/60, speed_y./sum(speed_y)); xlim([0 50]); hold on; 
% plot(speed_x*1000/60/60, speed_y./sum(speed_y)); hold on; 




end