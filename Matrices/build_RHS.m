%function build_RHS(FEmatrices,param)

RHSderiv = cell(1,5);

Xref = [0;0];

%background_medium = FEmatrices.Nodes(FEmatrices.background_nodes,:);

c0 = 340;

f = sym('f');
theta = sym('theta');
x = sym('x');
y = sym('y');

wave = @(f,theta,x,y) sin(2*pi*f/c0*((x-Xref(1))*cos(theta)+(y-Xref(2))*sin(theta)));


points = rand(2,100);

result = @(f,theta) 0;

tic;
for ii=1:length(points(1,:))
    result = @(f,theta) (result(f,theta) + wave(f,theta,points(1,ii),points(2,ii)));
end
toc;

%RHS(FEmatrices.background_nodes) = sin(2*pi*f);


%RHS = sparse(FEmatrices.size_system,1);


%end