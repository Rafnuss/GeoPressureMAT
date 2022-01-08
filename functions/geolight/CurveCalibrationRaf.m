function CurveCalibrationRaf(light,calib,known_coord)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


id = calib.first_period(1)<light.date&light.date<calib.first_period(2) | calib.second_period(1)<light.date&light.date<calib.second_period(2);

z = zenith(light.date, known_coord);

u = sqrt(6378/2/6.9)*sind(90-z);

c = 1; % cloudiness

SL = c*exp(-u.^2)./erfc(u);

SL(u<=0) = nan;
SL = SL/max(SL)*max(double(light.obs));


figure; hold on
plot(light.date,SL)
plot(light.date,light.obs)


figure; hold on
plot(light.date(id),log(SL))
plot(light.date(id),log(double(light.obs(id))))

figure; hold on
x=log(double(light.obs(id)));
x=movmedian(x,5);
x(isinf(x))=nan;
x(max(x)==x)=nan;
x(x<5)=nan;
y=log(SL(id));
y(hour(light.date(id))>5&hour(light.date(id))<15)=nan;


plot(y,x,'k.')

figure; hold on
plot(light.date(id),x)
plot(light.date(id),y)

plot(log(SL(id)),double(light.obs(id)),'k.')



end

