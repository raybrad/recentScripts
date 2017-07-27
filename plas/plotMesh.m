function plotMesh()

global nodes surfNodes;
figure;
patch('Vertices',nodes,'Faces',surfNodes,...
    'FaceVertexCData',hsv(1),'FaceColor','none')
view(-20, 30);