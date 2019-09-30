function stackedG = RealGFunction(state)
%G The observation function, takes a stacked state vector X and returns the
%observations for each point in the stacked state. n must equal 12
%(n elements per state).

%For the target object
% %Marker locations on object, [x,y,z; x2,y2,z2] in order
% %lb,blue,green,red,yellow
% Po = [0.125, 0.000, -0.070;
%     0.200, 0.100, -0.010;
%     0.000, 0.100, -0.010;
%     0.000, 0.000, -0.010;
%     0.075, 0.100, -0.070];

% %For the checkerboard target
% boardSize = [7,10];
% squareSize = 21.38;
% Po = generateCheckerboardPoints(boardSize,squareSize);
% Po = [Po, zeros(54,1)];
% Po = Po./1000;

%For the downsampled checkerboard target. Select points 1,6,28 (middle),49,54 
Po = [0,0,0; 0,106.9,0; 85.52,64.14,0; 171.04,0,0; 171.04,106.9,0];
Po = Po./1000;

% %For realsense
% camMat = [901.3116,         0,  645.9995;...
%     0,  905.3527,  358.4041;...
%     0,         0,    1.0000];

%For grasshopper
camMat = [852.73,         0,  614.05;...
    0,    857.10,    530.80;...
    0,         0,    1.0000];

sz = size(state,1);
assert(mod(sz,12) == 0); %State must have N*12 elements
sz = sz/12;

% if(size(state,1) ~= 12 || size(state,2) ~= 1)
%     disp(state);
%     disp("ERROR: bad state vector");
% end

%Function code ------------------------------------------
stackedG = [];
for step = 0:sz -1
    x = state(1 + step*12);
    y = state(3 + step*12);
    z = state(5 + step*12);
    a1 = state(7 + step*12);
    a2 = state(9 + step*12);
    a3 = state(11 + step*12);

    G_of_x = [];
    for n = 1:size(Po,1)

        %Put markers into camera frame in XYZ order
        R = [cos(a2)*cos(a3), -cos(a2)*sin(a3), sin(a2);
            cos(a1)*sin(a3)+cos(a3)*sin(a1)*sin(a2), cos(a1)*cos(a3)-sin(a1)*sin(a2)*sin(a3), -cos(a2)*sin(a1);
            sin(a1)*sin(a3)-cos(a1)*cos(a3)*sin(a2), cos(a3)*sin(a1)+cos(a1)*sin(a2)*sin(a3), cos(a1)*cos(a2)];

        T1 = [x;y;z];
        T2 = Po(n,:)';
        position = T1 + R*T2;

        %Does not work because compute_transform cannot take imaginary args
        %     tf = compute_transform(x,y,z,a1,a2,a3);
        %     point = [Po(n,:),1]';
        %     position = tf*point;


        Xc(n) = position(1);
        Yc(n) = position(2);
        Zc(n) = position(3);

        projection = camMat*[Xc(n)/Zc(n);Yc(n)/Zc(n);1];
        G_of_x = [G_of_x;projection(1:2)]; %Stack them

    end
    stackedG = [stackedG; G_of_x];
    
end

% disp(G_of_x);

% if(min(G_of_x) < 0)
%     disp("ERROR: NEGATIVE PIXEL VALUE FOUND");
%     disp(noVar);
% end

end
