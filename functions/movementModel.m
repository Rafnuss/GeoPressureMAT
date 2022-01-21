function mvt_pdf = movementModel(method, varargin)

switch method
    case "gam"
        if nargin>2
            mvt_pdf =  @(x) gampdf(x,varargin{1},varargin{2});
        else
            mvt_pdf =  @(x) gampdf(x,7,7);
        end
        
    case "norm"
        mvt_pdf = @(x) normpdf(c,30, 10);

    case "energy"
        bird = Bird(varargin{1});
        Pmech = mechanicalPower(bird);
        
        mvt_pdf_ms = matlabFunction((1./Pmech)^3/vpaintegral(1./Pmech,0,1000));
        mvt_pdf = @(x) mvt_pdf_ms(max(x,5)*1000/60/60);
    case "step"
        step = 50; rate=5;
        if nargin>2
            step = varargin{1};
        end
        if nargin>3
            rate = varargin{2};
        end
        mvt_pdf = @(x) min(exp(-(x-step)/rate),1);
end


% figure; 
if false
    speed_x=1:1000;
    figure;
    plot(speed_x, log(mvt_pdf(speed_x))); xlim([0 100]);
    %figure; plot(speed_x*1000/60/60, (mvt_pdf(speed_x))); xlim([0 100]);
end




end