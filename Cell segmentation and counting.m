%% Análise de Imagem Biomédica (EBE0056)
% Grupo 14 
% Ana Maria Sousa  up201707312
% Bernardo Pereira up201909521
% Jacinta Ferreira up201705451
%% Task2 segmentação e deteção de celulas 
close all, clear all;

%leitura de imagens, Roi e location de um diretório
imagefiles = dir ('*.tiff');
nfiles = length(imagefiles);
maskfiles = dir ('*.png');
locfiles = dir ('*.mat');

for i=1:nfiles
    
    %abrir imagen
    imagefilename = imagefiles(i).name;
    img0 = imread(imagefilename);
    
    % abrir ground truth
    maskfilename = maskfiles(i).name;
    mask = imread(maskfilename);
    
    % abrir localizações
    locfilename = locfiles(i).name;
    locs=load(locfilename);
    
    %conversão da imagem para tons de cinzento e double
    img1=rgb2gray(img0);
    img1=im2double(img1);
    %figure; imshow(img1); title('imagem inicial');
    
    %leitura da Roi
    %dilatação da ROI para facilitar deteção de celulas nos
    %limites esquerdo e superior
    se=strel('rectangle',[58,58]);  
    mask=imdilate(mask,se);
    img0=img1.*mask;
    %figure;imshow(img1); title('Imagem na ROi dilatada');
    
    %% Equalização adaptativa dos histograma (aumenta contraste)
    img1=adapthisteq(img0);
    
    %Filtro gaussiano 
    %(mascara gausiana de dimensão 4sigma+1)
    H = fspecial('gaussian',[5 5],1);    
    h2=imfilter(img1,H);
    %figure; imshow(h2); title('Filtro gaussiano');
    
    %determinação do gradiente
    [Gmag,Gdir] = imgradient(h2);
    %figure; imshow(Gmag, []); title('Magnitude');
    
    %% obter linhas brancas e celulas sobrepostas pretas através binarização 
    T = multithresh(img1, 6);
    img3=im2bw(img1,T(6));
    %figure; imshow(img3); title('linhas brancas-celulas pretas');
    
    %obter linhas brancas por binarização 
    T = multithresh(img0, 3);
    img4=im2bw(img0,T(3));
    %figure; imshow(img4); title('esqueleto da imagem');
    
    %% Limites esquerdo e superior 
    %ligeira open-close das linhas verticais e horizontais destacando celulas
    N=3;
    se_v=strel('line',N,0);
    se_h=strel('line',N,90);
    img_sep_new1=imopen(imclose(img3,se_v),se_v);
    img_sep_new2=imopen(imclose(img3,se_h), se_h);
    
    %limite esquerdo e superior
    img_sepnew=(img_sep_new1&img_sep_new2);
    %figure,imshow(img_sepnew); title('limite superior e inferior');
    
    %% Obter limites do quadrado ROI 
    N=51;
    se_v=strel('line',N,0);
    se_h=strel('line',N,90);
    img_sep_v=imopen(imclose(img4,se_v),se_v);
    img_sep_h=imopen(imclose(img4,se_h), se_h);
    %figure,imshow(img_sep_v);
    %figure,imshow(img_sep_h);
    
    %% alargar e eliminar possiveis desconectividades (devido a celulas)
    %dos limites da ROI
    s=33; 
    se=strel('square',s);
    img_sep_v=imdilate(img_sep_v,se);
    img_sep_h=imdilate(img_sep_h,se);
    
    %Quadrado ROI e interseções
    img_sep=img_sep_v&img_sep_h;
    %figure,imshow(img_sep);title('Quadrado ROI e interseções');
    %%
    %Para encontrar as limites
    img_sep1=bwhitmiss(img_sep,ones(25,100));
    img_sep2= bwhitmiss(img_sep, ones(100,25));
    img_sep3=img_sep1;
    img_sep3(1:600,1:1600)=0;
    img_sep1(600:1200,1:1600)=0;
    
    img_sep4=img_sep2;
    img_sep4(1:1200,1:800)=0;
    img_sep2(1:1200,800:1600)=0;
    
    %% binarização da magnitude do gradiente
    G=imbinarize(Gmag, 0.35); 
    %figure; imshow(G); title('Gmag após binarização');
    
    %% limite esquerdo e superior (celulas a branco, linhas a preto)
    img_sepnew=~img_sepnew.*(img_sep2|img_sep1);
    %figure; imshow(img_sepnew); title('limite superior e esquerdo');
    
    se1=strel('rectangle', [5 5]); 
    img_sepnew=imclose(imerode((img_sepnew),se1),se1);
    
    %% Construção da imagem da imagem combinada
    
    %limite supeior e esquerdo
    sep=(img_sep1|img_sep2);
    %figure; imshow(sep);
    
    %Construção da imagem combinada (Gradiente binarizado+limites 
    %esquerdo/superior após processamento)
    G=(G.*~(sep))|(img_sepnew);
    %figure; imshow(G); title('Imagem combinada');
    
    %% Dilatação dos limites direito e inferior ajustando a largura de 
    %forma a impedir a captação de celulas que toquem as linhas do meio
       se=strel('rectangle',[31 5]);  
       img_sep3=imdilate(img_sep3,se);
       se=strel('rectangle',[5 29]);
       img_sep4=imdilate(img_sep4,se);
    
    %% Preenchimento das linhas internas
    img3=img3-(img3&img_sep1)-(img3&img_sep2);
    %figure;imshow(img3,[]);
    
    G=G|img3;
    %figure;imshow(G,[]);
    
    %% Deteção de celulas
    [centers, radii, metrics]=imfindcircles((G|img_sep3|img_sep4),[23 51],'Sensitivity',0.885); 
    %figure;imshow(img0,[]);
    %viscircles(centers, radii,'EdgeColor','y');

    %% Avaliação
    
    tp(i)=0;
    cells=size(locs.cell_roi_pos);
    cellNumber(i)=cells(1);
    for j=1:cellNumber(i)
        location{j}=zeros(1200,1600);
        for a= locs.cell_roi_pos(j,2): locs.cell_roi_pos(j,2)+ locs.cell_roi_pos(j,4)
            for  b= locs.cell_roi_pos(j,1): locs.cell_roi_pos(j,1)+ locs.cell_roi_pos(j,3)
                location{j}(a,b)=1;
            end
        end
    end
    
    
    centers=round(centers);
    radii=floor(radii);
    foundNumber(i)=length(centers);
    for j=1:foundNumber(i)
        foundLocation{j}=zeros(1200,1600);
        for k=(centers(j,2)-radii(j)):(centers(j,2)+radii(j))
            for  l=(centers(j,1)-radii(j)):(centers(j,1)+radii(j))
                foundLocation{j}(k,l)=1;
            end
        end
        
    end
    for a=1:foundNumber(i)
        foundLocation{a}=foundLocation{a}(1:1200,1:1600);
        for b=1:cellNumber(i)
            location{b}=location{b}(1:1200,1:1600);
            intersection=bitand(foundLocation{a},location{b});
            jacc=sum(sum(intersection))/(sum(sum(foundLocation{a}))+sum(sum(location{b}))-sum(sum(intersection)));
            if jacc>=0.5
                tp(i)=tp(i)+1;
                break
            end
        end
    end
    fp(i)=foundNumber(i)-tp(i);
    fn(i)=cellNumber(i)-tp(i);
    r(i)=tp(i)/cellNumber(i);
    p(i)=tp(i)/foundNumber(i);
    f1(i)=2*p(i)*r(i)/(p(i)+r(i));
end

%calculo da media dos diversos criterios avaliados para as imagens
FoundNumber=mean(foundNumber)
TP=mean(tp)
FP=mean(fp)
FN=mean(fn)
R=mean(r)
P=mean(p)
F1=mean(f1)
