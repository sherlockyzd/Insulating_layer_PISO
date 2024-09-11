function [ s ] = Writer( u,v,p,xp,yp,data,Lx,Ly,XL,YL,NL,Nx,Ny,Iter )

% dx = data.dx;
dy = data.dy;
% hx = data.hx;
% hy = data.hy;
% vx = data.vx;
% vy = data.vy;
U_lid = data.U_lid;
% rho = data.rho;
% mu = data.mu;

s = zeros(Nx+2,Ny+2);
for i = 2:Nx+1
    for j = 2:Ny+1
        s(i,j) = s(i,j-1) + 0.5*(u(i,j)+u(i-1,j))*dy;
    end
end

Export = zeros((Nx+2)*(Ny+2),7);
k = 1;
i = 1;
for j = 1:Ny+2
    Export(k,1) = xp(i);
    Export(k,2) = yp(j);
    Export(k,3) = 0.0;
    Export(k,4) = 0.0;
    Export(k,5) = p(i,j);
    Export(k,6) = s(i,j);
    Export(k,7) = 0.0;
    k = k+1;
end
for i = 2:Nx+1
    j = 1;
    Export(k,1) = xp(i);
    Export(k,2) = yp(j);
    Export(k,3) = 0.0;
    Export(k,4) = 0.0;
    Export(k,5) = p(i,j);
    Export(k,6) = s(i,j);
    Export(k,7) = 0.0;
    k = k+1;
    for j = 2:Ny+1
        Export(k,1) = xp(i);
        Export(k,2) = yp(j);
        Export(k,3) = 0.5*(u(i,j)+u(i-1,j));
        Export(k,4) = 0.5*(v(i,j)+v(i,j-1));
        Export(k,5) = p(i,j);
        Export(k,6) = s(i,j);
        Export(k,7) = sqrt((0.5*(u(i,j)+u(i-1,j)))^2.0+(0.5*(u(i,j)+u(i-1,j)))^2.0);
        k = k+1;
    end
    j = Ny+2;
    Export(k,1) = xp(i);
    Export(k,2) = yp(j);
    Export(k,3) = U_lid;
    Export(k,4) = 0.0;
    Export(k,5) = p(i,j);
    Export(k,6) = s(i,j);
    Export(k,7) = U_lid;
    k = k+1;
end
i = Nx+2;
for j = 1:Ny+2
    Export(k,1) = xp(i);
    Export(k,2) = yp(j);
    Export(k,3) = 0.0;
    Export(k,4) = 0.0;
    Export(k,5) = p(i,j);
    Export(k,6) = s(i,j);
    Export(k,7) = 0.0;
    k = k+1;
end
fileID = fopen(['res1_',num2str(Iter),'.plt'],'w');
fprintf(fileID,'%s \n','TITLE = "PipeFlow1"');
fprintf(fileID,'%s \n','VARIABLES = "X" , "Y" , "U" , "V" , "P" , "S" , "Vm"');
fprintf(fileID,'%s %d %s %d %s %d %s \n',' ZONE I = ',Ny+2,', J = ',Nx+2,', K = ',1,', DATAPACKING=POINT');
fprintf(fileID,'%g %g %g %g %g %g %g\r\n',Export');
fclose(fileID);

X = zeros(Nx+2,Ny+2);
Y = zeros(Nx+2,Ny+2);
U = zeros(Nx+2,Ny+2);
V = zeros(Nx+2,Ny+2);
k = 1;
for i = 1:Nx+2
    for j = 1:Ny+2
        X(i,j) = Export(k,1);
        Y(i,j) = Export(k,2);
        U(i,j) = Export(k,3);
        V(i,j) = Export(k,4);
        k = k+1;
    end
end

Xq = XL*ones(1,NL+1);
Yq = 0.0:(Ly/NL):Ly;
Qq = interp2(X',Y',U',Xq,Yq);

fileID = fopen(['res2_',num2str(Iter),'.plt'],'w');
fprintf(fileID,'%s \n','TITLE = "PipeFlow2"');
fprintf(fileID,'%s \n','VARIABLES = "U" , "Y"');
for i = 1:NL+1
    fprintf(fileID,'%g %g \r\n',Qq(i),Yq(i));
end
fclose(fileID);

Xq = 0.0:(Lx/NL):Lx;
Yq = YL*ones(1,NL+1);
Qq = interp2(X',Y',V',Xq,Yq);

fileID = fopen(['res3_',num2str(Iter),'.plt'],'w');
fprintf(fileID,'%s \n','TITLE = "PipeFlow3"');
fprintf(fileID,'%s \n','VARIABLES = "X" , "V"');
for i = 1:NL+1
    fprintf(fileID,'%g %g \r\n',Xq(i),Qq(i));
end
fclose(fileID);


end

