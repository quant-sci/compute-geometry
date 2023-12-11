%------------------------------ PolyScript -------------------------------%
% Chama a Fun��o PolyMesher com os seguintes par�metros:                  %
% @MBBDomain = Tipo de Dom�nio = VIGA MBB                                 %
% 5000       = N�mero total de Pol�gonos a serem gerados                  %
% 50         = N�mero de Itera��es de Lloyd (para fazer com que os        %
%              pol�gonos fiquem quase "regulares")
%-------------------------------------------------------------------------%

%% ---------------------------------------------------- CREATE 'fem' STRUCT
[Node,Element,Supp,Load] = PolyMesher(@MbbDomain,1000,100);

%% ------------------------------------------------------------------------