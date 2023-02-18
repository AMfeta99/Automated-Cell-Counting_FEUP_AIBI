%% Análise de Imagem Biomédica (EBE0056)
% Grupo 14
% Ana Maria Sousa  up201707312
% Bernardo Pereira up201909521
% Jacinta Ferreira up201705451
%% Task 1- Determinação da região de interesse (ROI)

%abrir imagens
imagefiles = dir ('*.tiff');
nfiles = length(imagefiles);
for i=1:nfiles
    currentfilename = imagefiles(i).name;
    currentimage = imread(currentfilename);
    images{i} = currentimage;
end

% abrir ground truth
maskfiles = dir ('*.png');
for i=1:nfiles
    currentfilename = maskfiles(i).name;
    currentmask = imread(currentfilename);
    masks{i} = currentmask;
end

for i=1:nfiles
    
    %Leitura da imagem
    img= images{i};
    mask = masks{i};
    img=rgb2gray(img); %converter de RGB para tons de cinzento
    %figure,imshow(img1, []); title('imagem inicial');
    
    %Colocar imagem a preto e branco
    T = graythresh(img);
    imgBin2=im2bw(img,T);
    %figure,imshow(imgBin2,[]); title('imagem binarizada usando threshold');
    
    %% Percorrer imagem
    % achar linhas
    %vetores que representa as arestas horizontais do quadrado
    X=ones(1600,1);
    X2=ones(1600,1);
    
    for jj=1:1600
        linha= false;
        linha2= false;
        
        %Percorrer a primeira metade da imagem ou até achar ponto que pertence à linha
        for ii=1:600
            %primeira linha
            if (imgBin2(ii,jj)==1 && linha==false)
                
                if (imgBin2(ii+1,jj)==0)
                    k=ii+15;
                    
                    if (imgBin2(k,jj)==1)
                        x=ii-10;
                        linha=true;
                        X(jj,1)=x;
                    end
                end
            end
        end
        
        %Percorrer a segunda metade da imagem ou até achar ponto que pertence à linha
        for ii=1200:-1:600
            %ultima linha
            if (imgBin2(ii,jj)==1 && linha2==false)
                
                if (imgBin2(ii-1,jj)==0)
                    k2=ii-15;
                    
                    if (imgBin2(k2,jj)==1)
                        x2=ii+10;
                        linha2=true;
                        X2(jj,1)=x2;
                    end
                end
            end
        end
    end
    
    %ordenar vetores
    X=sort(X);
    X2=sort(X2);
    
    %% achar colunas
    %vetores que representa as arestas verticais do quadrado
    Y=ones(1200,1);
    Y2=ones(1200,1);
    
    %Percorrer a primeira metade da imagem ou até achar ponto que pertence à
    %ROI
    for ii=1:1200
        coluna= false;
        coluna2= false;
        
        for jj=1:800
            %primeira coluna
            if (imgBin2(ii,jj)==1 && coluna==false) %encontrar as 3 colunas juntas
                
                if (imgBin2(ii,jj+1)==0)
                    k=jj+15;
                    
                    if (imgBin2(ii,k)==1)
                        y=jj-10;                        %dimensão da barra é ~10
                        coluna=true;
                        Y(ii,1)=y;
                    end
                end
            end
        end
        
        %Percorrer a primeira metade da imagem ou até achar ponto que pertence à
        %aresta
        for jj=1600:-1:800
            %Ultima coluna
            if (imgBin2(ii,jj)==1 && coluna2==false)
                
                if (imgBin2(ii,jj-1)==0)
                    k2=jj-15;
                    
                    if (imgBin2(ii,k2)==1)
                        y2=jj+10;
                        coluna2=true;
                        Y2(ii,1)=y2;
                        
                    end
                end
            end
        end
        
    end
    %ordenar vetores
    Y=sort(Y);
    Y2=sort(Y2);

imgBin3=imgBin2;
%recortar imagem assumindo com retas a que passa pela mediana 
imgBin3(1:median(X),1:1600)=0;
imgBin3(median(X2):1200,1:1600)=0;
%figure,imshow(imgBin2); title('imagem parte superior/inferior');

%usar mediana como valor para as retas 
imgBin2(1:1200,1:median(Y))=0;
imgBin2(1:1200,median(Y2):1600)=0;

%interseção das imagens
imgBin2=imgBin2&imgBin3;
%figure, imshow(imgBin2); title('imagem Roi');   

% quadrado final
imgBin2(median(X):median(X2),median(Y):median(Y2))=1;
% figure; imshow(imgBin2); title('ROI');
    
%% Indice de jaccard
   
 intersection=bitand(imgBin2, mask);
 jacc(i)=sum(sum(intersection))/(sum(sum(imgBin2))+sum(sum(mask))-sum(sum(intersection)));
    
%% Distância entre os vertices 
%encontrar os vertices
[I,J]=find(mask>max(mask(:))/2);
IJ=[I,J];
[~,idx]=min(IJ*[1 1; -1 -1; 1 -1; -1 1].');
corners=IJ(idx,:);
    
% ordenar os vertices
for b=1:4
    if corners(b,1)<600 && corners(b,2)<800
        ordCorners(1,:)=corners(b,:);
    elseif corners(b,1)<600 && corners(b,2)>800
        ordCorners(2,:)=corners(b,:);
    elseif corners(b,1)>600 && corners(b,2)<800
        ordCorners(3,:)=corners(b,:);
    elseif corners(b,1)>600 && corners(b,2)>800
        ordCorners(4,:)=corners(b,:);
    end
end
    
%encontrar as distancias
dist(i,1)=sqrt(((ordCorners(1,1)-median(X))^2)+ ((ordCorners(1,2)-median(Y))^2));
dist(i,2)=sqrt(((ordCorners(2,1)-median(X))^2)+ ((ordCorners(2,2)-median(Y2))^2));
dist(i,3)=sqrt(((ordCorners(3,1)-median(X2))^2)+ ((ordCorners(3,2)-median(Y))^2));
dist(i,4)=sqrt(((ordCorners(4,1)-median(X2))^2)+ ((ordCorners(4,2)-median(Y2))^2));

%média e máximo das distâncias
maxDist(i)=max(dist(i,:));
meanDist(i)=mean(dist(i,:));
end

%% Médias dos critérios de avaliação
jaccFinal = mean (jacc);
meanDistFinal = mean(meanDist);
maxDistFinal= mean(maxDist);



