function [repro_error, x2_repro] = backprojection(correspondences, P1, Image2, T, R, K)
    % This function calculates the mean error of the back projection
    % of the world coordinates P1 from image 1 in camera frame 2
    % and visualizes the correct feature coordinates as well as the back projected ones.
    dim = size(P1);
    
    % Weltkoordinaten der Punkte P1 bezogen auf Kamerasystem 2
    P2 = R*P1 + repmat(T,1,dim(2));
    
    % Projektion der Punkte P2 auf die Bildebene 2
    P2_proj = P2./P2(3,:);
    
    % Umwandlung der Bildkoordinaten in Pixelkoordinaten
    x2_tilde = K*P2_proj;
    
    % Image2 anzeigen lassen
    figure;
    hold on
    imshow(Image2);
    
    % x2_tilde und x2 plotten und beschriften
    plot(x2_tilde(1,:),x2_tilde(2,:),'xgreen');
    plot(correspondences(3,:),correspondences(4,:),'xred');
    text(x2_tilde(1,:),x2_tilde(2,:),string([1:dim(2)]),'Color','green','FontSize',8);
    text(correspondences(3,:),correspondences(4,:),string([1:dim(2)]),'Color','red','FontSize',8);
    
    % Berechnung des durchschnittlichen Fehlers
    error = 0;
    for i = 1:dim(2)
        error = error + norm(correspondences(3:4,i) - x2_tilde(1:2,i));
    end
    mean_error = error/dim(2);
    
    % Variablenzuweisung
    repro_error = mean_error;
    x2_repro = x2_tilde;
end