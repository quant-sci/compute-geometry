%------------------------------ PolyScript -------------------------------%
% Chama a Funcao PolyMesher com os seguintes parametros:                  %
% @MBBDomain = Tipo de Dom�nio = VIGA MBB                                 %
% 1500       = N�mero total de Pol�gonos a serem gerados                  %
% 100         = N�mero de Itera��es de Lloyd (para fazer com que os        %
%              pol�gonos fiquem quase "regulares")
%-------------------------------------------------------------------------%

%% ---------------------------------------------------- CREATE 'fem' STRUCT
[Node,Element,Supp,Load] = PolyMesher(@MbbDomain,1500,100,[2.5 0.1],500,[0.1 1],1,[0.8 1 0 0]);

%% ------------------------------------------------------------------------