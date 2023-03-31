
%% TRANSATMOSPHERIC ASSIGNMENT %%%%
%% Carleton University %%%%
%% Tyler Blair - 101074740 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all


dTdh = [-6.5 0 1 2.8 0 -2.8 -2]; 
atmLevels  = [0 11 20 32 47 51 72 86]; %km

ag = 9.81; %m/s^2
R = 287; %J/kgK
K = 1.4;

atmDensity(1) = 1.225; %Kg/m^3
T(1) = 288.15; % K
a(1) = sqrt(K*R*T(1));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Creates a vector of atmospheric densities and temps at each altitude level
for i = 1:length(dTdh)

    value = ag/(R*dTdh(i)/1000);
    levelDiff = (atmLevels(i+1)-atmLevels(i))*1000;

    if dTdh(i) ~= 0
        atmDensity(i+1) = atmDensity(i)*((1+(((dTdh(i))*((atmLevels(i+1)-atmLevels(i))))/(T(i))))^(-1*(value+1)));
    else
        atmDensity(i+1) = atmDensity(i)*(exp(-1*(ag*(levelDiff)/(R*T(i)))));
    end

  T(i+1) = T(i)+dTdh(i)*(levelDiff/1000);
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rocket %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

launchAngles = (15:2.5:25);
payloads = (225:100:525);

for i = 1:length(launchAngles)
    for k = 1:length(payloads)
        
        theta = deg2rad(launchAngles(i));
        mPayload = payloads(k);

        isp = 225; %s
        mProp = 1015; % Kg
        mMotor = 252; % Kg
        mStructure = mPayload*0.5;
        totalMass = mProp+mMotor+mPayload+mStructure;
        rocketLength = 11.38; %m
        rocketDiameter = 0.4384; %m

        gamma = totalMass / (totalMass - mPayload);
        beta = mStructure / (totalMass - mPayload);
        mu = (gamma + 1) / (gamma + beta);

        

        

        rE = 6378000; % Earth radius km
        sRef = (pi()*rocketDiameter^2)/4;
        sWet = pi()*rocketDiameter*rocketLength;

        %% Rocket Thrust %%
        rocketTime  = [0 0.22 0.33 0.77 2.16 3.15 5.13 8.10 12.06 16.02 19.97 21.95 24.92 26.90 28.19 28.88 30.07 31.11 31.23 5000];
        rocketThrust = [0 79419 77131 74838 73353 73413 75105 75287 78662 78914 76014 76136 73141 70875 68538 61354 46156 1805 0 0] ;

        %% Sets up time vector for simulation %%
        timeStep = 0.05;
        time = linspace(0,rocketTime(end),(rocketTime(end)/timeStep)+1);



        %% Rocket Trajectory %%
        vRocket = 0;
        mach = 0;
        rRocket = 6378; %Rocket starting total altitude km


        deltaV(1) = 0;
        instMass(1) = totalMass;
        alt(1) = 0;
        range(1) = 0;
        phi(1) = 0;
        Machspeed(1) = 0;

        thrust = zeros(1,length(time));
        vRocket = zeros(1,length(time));
        dynP = zeros(1,length(time));
        pDrag = zeros(1,length(time));
        frDrag = zeros(1,length(time));


        isBurn = 1;
        isBallistic = 0;
        inAtmosphere = 1;
        j = 1;
        totalPropMass = mProp;



        %% Rocket Burn phase %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        while isBurn == 1


        %% Linear interpolation for thrust %%
        if time(j) < rocketTime(end) && totalPropMass >= 0

            thrustIndex = find(rocketTime > time(j), 1) - 1;
            thrust(j) = (time(j) - rocketTime(thrustIndex)) * ((rocketThrust(thrustIndex + 1) - rocketThrust(thrustIndex)) / (rocketTime(thrustIndex + 1) - rocketTime(thrustIndex))) + rocketThrust(thrustIndex);

        else 

            thrust(j) = 0;

        end

            %% Calculates current propellant massflow and total mass %%
            massflowProp = thrust(j)/(isp * ag); 
            totalPropMass = totalPropMass - massflowProp * (time(j+1) - time(j));
            instMass(j+1) = instMass(j) - massflowProp * (time(j+1) - time(j)); 

            %% Determines if rocket is still in burn phase %%
        if totalPropMass <= 0
            isBurn = 0;
            burnTime = time(j);
            burnIndex = j;
            isBallistic = 1;
            ballisticMass = instMass(j+1);
            break
        end


            %% Gets parameters associated with atmospheric model %%
            [density,mach,temp,cD,cFr] = atmoCalc((alt(j)/1E3),vRocket(j),T,atmDensity);
            
            %% Calculates drag forces %%
            dynP(j) = 0.5*density*vRocket(j)^2; % Dynamic pressure
            pDrag(j) = cD*dynP(j)*sRef;         % Pressure drag
            frDrag(j) = (cFr*dynP(j)*sRef);     % Friction drag


            %% Calculates change in velocity %%
            
                deltaV(j) = ((isp*ag*(instMass(j)-instMass(j+1)))/instMass(j+1)) + (-(pDrag(j)/instMass(j+1))-(frDrag(j)/instMass(j+1))-ag*cos(theta))*(time(j+1)-time(j)); 
            
        

            %% Calculates new velocity and altitude %%

            if alt(j) == 0 && deltaV(j) <= 0
            vRocket(j+1) = 0;
            else
            vRocket(j+1) = vRocket(j)+deltaV(j);
            end
            Machspeed(j+1) = mach;

            %% Calculates new rocket position %%
            distanceRad = vRocket(j+1)*(time(j+1)-time(j))*cos(theta)+rE; % radial distance
            distanceTan = vRocket(j+1)*(time(j+1)-time(j))*sin(theta);    % tangential distance
            distanceHyp = sqrt(distanceRad^2+distanceTan^2);              % hypotenuse distance
            phi(j+1) = phi(j)+atan(distanceTan/distanceRad);                       % Angular distance of rocket from center of earth
            range(j+1) = phi(j)*rE;                            % Circular distance of rocket from center of earth
            alt(j+1) = alt(j)+(distanceRad-rE);
        


            j = j+1;
            
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %% End of burn phase %%

        


        %% Rocket Ballistic phase %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %% Sets x and y components of velocity %%
        vRocketX(j) = vRocket(j)*sin(theta);
        vRocketY(j) = vRocket(j)*cos(theta);


        while isBallistic == 1 || time(j) > 1000

            if inAtmosphere == 1
        
                %% Gets parameters associated with atmospheric model %%
                [density,mach,temp,cD,cFr] = atmoCalc((alt(j)/1E3),vRocket(j),T,atmDensity);
            
                %% Calculates drag forces %%
                dynP(j) = 0.5*density*vRocket(j)^2; % Dynamic pressure
                pDrag(j) = cD*dynP(j)*sRef;         % Pressure drag
                frDrag(j) = (cFr*dynP(j)*sRef);     % Friction drag

        
            else  
                %% If rocket is out of atmosphere, drag forces are zero %%
                dynP(j) = 0;
                pDrag(j) = 0; 
                frDrag(j) = 0;
            end

            %% Calculates change in velocity %%
            deltaV(j) = (-(pDrag(j)/ballisticMass)-(frDrag(j)/ballisticMass)-ag*cos(theta))*(time(j+1)-time(j));  
        
            %% Calculates new velocity %%
            vRocketX(j+1) = vRocketX(j)+(-((pDrag(j)/ballisticMass)-(frDrag(j)/ballisticMass))*(time(j+1)-time(j))*sin(theta));
            vRocketY(j+1) = vRocketY(j)+(-((pDrag(j)/ballisticMass)-(frDrag(j)/ballisticMass))*(time(j+1)-time(j))*cos(theta))-ag*(time(j+1)-time(j));
            vRocket(j+1) = sqrt((vRocketX(j+1))^2+(vRocketY(j+1))^2);
            Machspeed(j+1) = mach;

        
            %% Calculates new rocket position %%
            distanceRad = rE+vRocket(j+1)*(time(j+1)-time(j))*cos(theta); % radial distance
            distanceTan = vRocket(j+1)*(time(j+1)-time(j))*sin(theta);    % tangential distance
            distanceHyp = sqrt(distanceRad^2+distanceTan^2);              % hypotenuse distance
            phi(j+1) = phi(j)+atan(distanceTan/distanceRad);              % Angular distance of rocket from center of earth
            range(j+1) = phi(j+1)*rE;                               % Circular distance of rocket from center of earth
            alt(j+1) = alt(j)+distanceRad-rE;
            
        

            %% Calculates new rocket angle using projectile motion %%
            if vRocketY(j+1) > 0
                theta = atan(vRocketX(j+1)/vRocketY(j+1));
            elseif vRocketY(j+1) < 0
                theta = deg2rad(90)+atan(abs(vRocketY(j+1))/vRocketX(j+1));
            else
                theta = deg2rad(90);
            end
        




            %% Checks if rocket makes ground contact %%
            if alt(j+1) <= 0
                alt(j+1) = 0;
                isBallistic = 0;
                endindex = j;
                break
            end

            %% Checks if rocket is still in atmosphere %%
            if (alt(j+1)/1E3) < atmLevels(end)
                inAtmosphere = 1;
            else
                inAtmosphere = 0;
            end

            j = j+1; 

        end

        %% Caculates Apogee and time of Apogee %%
        apogee = max(alt);
        apogeeTime = time(find(alt == apogee));
        apogeeIndex = find(alt == apogee);


        %% Stores data for each payload %%

        payloadData(k).payload = mPayload;
        payloadData(k).time = time;
        payloadData(k).alt = alt;
        payloadData(k).range = range;
        payloadData(k).vRocket = vRocket;
        payloadData(k).vRocketX = vRocketX;
        payloadData(k).vRocketY = vRocketY;
        payloadData(k).rocketMass = instMass;
        payloadData(k).theta = theta;
        payloadData(k).phi = phi;
        payloadData(k).dynP = dynP;
        payloadData(k).pDrag = pDrag;
        payloadData(k).frDrag = frDrag;
        payloadData(k).deltaV = deltaV;
        payloadData(k).Machspeed = Machspeed;
        payloadData(k).burnTime = burnTime;
        payloadData(k).burnIndex = burnIndex;
        payloadData(k).apogee = apogee;
        payloadData(k).apogeeTime = apogeeTime;
        payloadData(k).apogeeIndex = apogeeIndex;
        payloadData(k).endIndex = endindex;

        payloadData(k).Simple.avgThrust = mean(rocketThrust(1:end-2))/9.81;
        
        payloadData(k).Simple.vBurnout(1) = log(mu)*(isp * ag);
        payloadData(k).Simple.timeBurnout = (isp/(payloadData(k).Simple.avgThrust/totalMass))*(1-(1/mu));
        payloadData(k).Simple.vBurnout(2) = isp*ag*log(mu)-ag*(cos(theta))*payloadData(k).Simple.timeBurnout;
        payloadData(k).Simple.hBurnout = (isp*ag)*payloadData(k).Simple.timeBurnout*(1-(log(mu)/(mu-1)))-0.5*ag*payloadData(k).Simple.timeBurnout^2;
        payloadData(k).Simple.hMax = payloadData(k).Simple.hBurnout + payloadData(k).Simple.vBurnout(2)^2/(2*ag);

       







    end

    launchAngleData(i).launchAngle = theta;
    launchAngleData(i).payloadData = payloadData;
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% End of ballistic phase %%

%% Creates menu asking to animate or plot %%

displayType = menu('Choose display type:','Animate','Plot');

if displayType == 1

    %% asks which payload and launch angle to animate from inputted arrays%%

    %% Display payload and launch angle arrays %%



    %% Asks user to input payload and launch angle to animate %%

    for i = 1:length(launchAngles)
        launchAnglesArray{i} = num2str(launchAngles(i));
    end

    for i = 1:length(payloads)
        payloadsArray{i} = num2str(payloads(i));
    end

    payloadAnimate = menu('Choose payload to animate:',payloadsArray);
    thetaAnimate = menu('Choose launch angle to animate:',launchAnglesArray);





    %% stores data for payload and launch angle to animate %%
    payloadData = launchAngleData(thetaAnimate).payloadData(payloadAnimate);
    time = payloadData.time;
    alt = payloadData.alt;
    range = payloadData.range;
    vRocket = payloadData.vRocket;
    instMass = payloadData.rocketMass;
    theta = payloadData.theta;
    phi = payloadData.phi;
    dynP = payloadData.dynP;
    pDrag = payloadData.pDrag;
    frDrag = payloadData.frDrag;
    deltaV = payloadData.deltaV;
    Machspeed = payloadData.Machspeed;
    burnTime = payloadData.burnTime;
    burnIndex = payloadData.burnIndex;
    apogee = payloadData.apogee;
    apogeeTime = payloadData.apogeeTime;
    apogeeIndex = payloadData.apogeeIndex;
    endIndex = payloadData.endIndex;


    %% animate rocket flight %%
    figure(2)
    for i = 1:8:endIndex
        
    
        %% plot trajectory
        if time(i) <= burnTime
        plot(range(1:i)/1000,alt(1:i)/1000,'r')
        hold on
        else
        plot(range(1:burnIndex)/1000,alt(1:burnIndex)/1000,'r')
        hold on
        plot(range(burnIndex:i)/1000,alt(burnIndex:i)/1000,'k')
        end
    
        %% plot rocket position
        plot(range(i)/1000,alt(i)/1000,'b.','MarkerSize',15)
    
        title('Rocket Trajectory')
        xlabel('Range (km)')
        ylabel('Altitude (km)')
        grid on


        %% fix axes to entire trajectory
        axis([-10 range(end)/1000+80 0 max(alt)/1000+10])

        %% adds text to right side to display all data, and shortens all values 

        text(range(end)/1000+20,1*max(alt)/1000,['Time: ' num2str(time(i),'%.1f') ' s'])
        text(range(end)/1000+20,0.9*max(alt)/1000,['Altitude: ' num2str(alt(i)/1000,'%.1f') ' km'])
        text(range(end)/1000+20,0.8*max(alt)/1000,['Range: ' num2str(range(i)/1000,'%.1f') ' km'])
        text(range(end)/1000+20,0.7*max(alt)/1000,['Velocity: ' num2str(vRocket(i),'%.1f') ' m/s'])
        text(range(end)/1000+20,0.6*max(alt)/1000,['Mach: ' num2str(Machspeed(i),'%.1f') ' m/s'])
        text(range(end)/1000+20,0.5*max(alt)/1000,['Drag: ' num2str(pDrag(i)/1000,'%.2f') ' kN'])
        text(range(end)/1000+20,0.4*max(alt)/1000,['Friction: ' num2str(frDrag(i)/1000,'%.2f') ' kN'])
        [density,~,~,cD,cFr] = atmoCalc((alt(i)/1E3),vRocket(i),T,atmDensity);
        text(range(end)/1000+20,0.3*max(alt)/1000,['Density: ' num2str(density,'%.2f') ' kg/m^3'])
        text(range(end)/1000+20,0.2*max(alt)/1000,['cD: ' num2str(cD,'%.2f')])
        text(range(end)/1000+20,0.1*max(alt)/1000,['cFr: ' num2str(cFr,'%.2f')])

        

            %% Label Burnout 
        if time(i) > burnTime
            plot(range(burnIndex)/1000,alt(burnIndex)/1000,'k.','MarkerSize',5)
            %% label altitude at burnout
            text(range(burnIndex)/1000,alt(burnIndex)/1000,['Burnout: ' num2str(alt(burnIndex)/1000,'%.1f') ' km'])
        end

            %% Label Apogee
    
        if time(i) > apogeeTime
            plot(range(apogeeIndex)/1000,alt(apogeeIndex)/1000,'k.','MarkerSize',5)
            %% label altitude at apogee
            text(range(apogeeIndex)/1000,alt(apogeeIndex)/1000,['Apogee: ' num2str(alt(apogeeIndex)/1000,'%.1f') ' km'])
        end
            %% label impact

        if alt(i) <= 0 && range(i) >= max(range)
            plot(range(i)/1000,alt(i)/1000,'k.','MarkerSize',5)
            text(range(i)/1000,alt(i)/1000,['Impact: ' num2str(alt(i)/1000,'%.1f') ' km'])
        end

        hold off
        pause(1/3600)
    end

    
elseif displayType == 2

    %% ask user to select between payloads, or for launch angles to createa composite plot
    plotType = menu('Select Plot Type','Payloads','Launch Angles');

    if plotType == 1

       
        %% ask for which launch angle to plot payload data
        %% launch angle array for menu
        for i = 1:length(launchAngles)
            launchAnglesString{i} = num2str(launchAngles(i));
        end   
        launchAngleSelect = menu('Select Launch Angle',launchAnglesString);

        for i = 1:length(payloads)
            payloadLegend{i} = ['Payload: ' num2str(payloads(i)) ' kg'];
        end

        
        

    
        for i = 1:length(payloads)
            payloadData = launchAngleData(launchAngleSelect).payloadData(i);
            endIndex = payloadData.endIndex;
            time = payloadData.time;
            alt = payloadData.alt;
            range = payloadData.range;
            vRocket = payloadData.vRocket;
   
            

            figure(1);

            %% plot altitude vs time
            subplot(3,1,1)
            plot(time(1:endIndex),alt(1:endIndex)/1000,'LineWidth',0.5)
            title('Altitude vs Time')
            xlabel('Time (s)')
            ylabel('Altitude (km)')
            grid on
            hold on

            %% plot range vs time
            subplot(3,1,2)
            plot(time(1:endIndex),range(1:endIndex)/1000,'LineWidth',0.5)
            title('Range vs Time')
            xlabel('Time (s)')
            ylabel('Range (km)')
            grid on
    
            hold on

            subplot(3,1,3)
            plot(time(1:endIndex),vRocket(1:endIndex),'LineWidth',0.5)
            title('Velocity vs Time')
            xlabel('Time (s)')
            ylabel('Velocity (m/s)')
            grid on
           hold on
            %% plot velocity vs time

            legend(payloadLegend)
            sgtitle(['Launch Angle: ' num2str(launchAngles(launchAngleSelect)) ' deg'])

        end

        for i = 1:length(payloads)

            payloadData = launchAngleData(launchAngleSelect).payloadData(i);
            endIndex = payloadData.endIndex;
            time = payloadData.time;
            pDrag = payloadData.pDrag;
            frDrag = payloadData.frDrag;
            dynP = payloadData.dynP;
            figure(2);
            %% plot drag , friction drag, and dynamic pressure vs time

            subplot(3,1,1)
            plot(time(1:endIndex),pDrag(1:endIndex)/1000,'LineWidth',0.5)

            title('Pressure Drag vs Time')
            xlabel('Time (s)')
            ylabel('Pressure Drag (kN)')
            grid on

            hold on

            subplot(3,1,2)
            plot(time(1:endIndex),frDrag(1:endIndex)/1000,'LineWidth',0.5)
            title('Friction Drag vs Time')
            xlabel('Time (s)')
            ylabel('Friction Drag (kN)')
            grid on

            hold on

            subplot(3,1,3)
            plot(time(1:endIndex),dynP(1:endIndex)/1000,'LineWidth',0.5)
            title('Dynamic Pressure vs Time')
            xlabel('Time (s)')
            ylabel('Dynamic Pressure (kPa)')
            grid on

            hold on

            legend(payloadLegend)
            sgtitle(['Launch Angle: ' num2str(launchAngles(launchAngleSelect)) ' deg'])


    
        end


    
        hold off


    elseif plotType == 2

    
        %% ask for which payload to plot launch angle data
        %% payload array for menu
        for i = 1:length(payloads)
            payloadsString{i} = num2str(payloads(i));
        end   
        payloadSelect = menu('Select Payload',payloadsString);

        for i = 1:length(launchAngles)
            launchAngleLegend{i} = ['Launch Angle: ' num2str(launchAngles(i)) ' deg'];
        end

        for i = 1:length(launchAngles)
            payloadData = launchAngleData(i).payloadData(payloadSelect);
            endIndex = payloadData.endIndex;
            time = payloadData.time;
            alt = payloadData.alt;
            range = payloadData.range;
            vRocket = payloadData.vRocket;
            pDrag = payloadData.pDrag;
            frDrag = payloadData.frDrag;
            dynP = payloadData.dynP;

            figure(1);
            %% plot altitude vs time
            subplot(3,1,1)
            plot(time(1:endIndex),alt(1:endIndex)/1000,'LineWidth',0.5)
            title('Altitude vs Time')
            xlabel('Time (s)')
            ylabel('Altitude (km)')
            grid on
            hold on

            %% plot range vs time
            subplot(3,1,2)
            plot(time(1:endIndex),range(1:endIndex)/1000,'LineWidth',0.5)
            title('Range vs Time')
            xlabel('Time (s)')
            ylabel('Range (km)')
            grid on
            hold on

            subplot(3,1,3)
            plot(time(1:endIndex),vRocket(1:endIndex),'LineWidth',0.5)
            title('Velocity vs Time')
            xlabel('Time (s)')
            ylabel('Velocity (m/s)')
            grid on
           hold on
            %% plot velocity vs time

            legend(launchAngleLegend)
            sgtitle(['Payload: ' num2str(payloads(payloadSelect)) ' kg'])


            figure(2);
            %% plot drag , friction drag, and dynamic pressure vs time

            subplot(3,1,1)
            plot(time(1:endIndex),pDrag(1:endIndex)/1000,'LineWidth',0.5)

            title('Pressure Drag vs Time')
            xlabel('Time (s)')
            ylabel('Pressure Drag (kN)')
            grid on
            hold on

            subplot(3,1,2)
            plot(time(1:endIndex),frDrag(1:endIndex)/1000,'LineWidth',0.5)
            title('Friction Drag vs Time')
            xlabel('Time (s)')
            ylabel('Friction Drag (kN)')
            grid on
            hold on

            subplot(3,1,3)
            plot(time(1:endIndex),dynP(1:endIndex)/1000,'LineWidth',0.5)
            title('Dynamic Pressure vs Time')
            xlabel('Time (s)')
            ylabel('Dynamic Pressure (kPa)')
            grid on
            hold on

            legend(launchAngleLegend)
            sgtitle(['Payload: ' num2str(payloads(payloadSelect)) ' kg'])
        end

    end

end

%% Plot density vs altitude

figure(3)
for i = 0:1:atmLevels(end)-1
    [density(i+1),mach,temp(i+1),cD,cFr] = atmoCalc(i,0,T,atmDensity);


end

plot(0:1:atmLevels(end)-1,density)
hold on

xlabel('Altitude (km)')
ylabel('Density (kg/m^3)')
grid on

yyaxis right
plot(0:1:atmLevels(end)-1,temp)

ylabel('Temperature (K)')
legend('Density','Temperature')
title('Density and Temperature vs Altitude')





function [density,mach,temp,cD,cFr] = atmoCalc(alt,vRocket,T,atmDensity)

    % This function calculates the density and speed of sound at a given altitude
    %% alt is the altitude in km

 
    dTdh = [-6.5 0 1 2.8 0 -2.8 -2]; 
    atmLevels  = [0 11 20 32 47 51 72 86];
    ag = 9.81;
    R = 287; %J/kgK
    K = 1.4;

   %% Atmosphere Check

  
   if alt == 0
       index = 1;
   elseif alt >= atmLevels(end)
     density = 0;
     mach = 0;
     temp = 0;
     cD = 0;
     cFr = 0;
     return
   else
    index = find(atmLevels>(alt),1)-1;
   end
    temp = T(index)+dTdh(index)*((alt)-atmLevels(index));

    value = ag/(R*(dTdh(index)/1000));

    if dTdh(index) ~= 0
        density = atmDensity(index)*((1+(((dTdh(index))*((alt-atmLevels(index))))/(T(index))))^(-1*(value+1)));
    else
        density = atmDensity(index)*(exp(-1*(ag*((alt-atmLevels(index))*1000)/(R*T(index)))));
    end

    
    vSound = sqrt(K*R*temp);
    mach = vRocket/vSound;

    if mach <= 0.5
        cD = 0.0798;
    elseif mach > 0.5 && mach <= 1.2
        cD = 0.394*(mach^(2.303));
    elseif mach > 1.2
        cD = 0.639*mach^(-0.488);
        
    end

    if mach >= 0.1 && mach <= 9
        cFr = (-0.0721)*log(mach)+0.1831;
    else
        cFr = 0;
    end

    


end