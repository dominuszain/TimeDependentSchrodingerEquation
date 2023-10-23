// This GUI file is generated by guibuilder version 4.2.2
//////////
f=figure('figure_position',[400,50],'figure_size',[1079,736],'auto_resize','on','background',[34],'figure_name','Graphic window number %d','color_map',[0,0,0;0,0,1;0,1,0;0,1,1;1,0,0;1,0,1;1,1,0;1,1,1;0,0,0.5647059;0,0,0.6901961;0,0,0.8156863;0.5294118,0.8078431,1;0,0.5647059,0;0,0.6901961,0;0,0.8156863,0;0,0.5647059,0.5647059;0,0.6901961,0.6901961;0,0.8156863,0.8156863;0.5647059,0,0;0.6901961,0,0;0.8156863,0,0;0.5647059,0,0.5647059;0.6901961,0,0.6901961;0.8156863,0,0.8156863;0.5019608,0.1882353,0;0.627451,0.2509804,0;0.7529412,0.3764706,0;1,0.5019608,0.5019608;1,0.627451,0.627451;1,0.7529412,0.7529412;1,0.8784314,0.8784314;1,0.8431373,0;0.8,0.8,0.8;0.9333333,0.9333333,0.9333333],'dockable','off','infobar_visible','off','toolbar_visible','off','menubar_visible','off','default_axes','on','visible','off');
//////////
handles.dummy = 0;
handles.tdse_axes= newaxes();handles.tdse_axes.margins = [ 0 0 0 0];handles.tdse_axes.axes_bounds = [0.0490104,0.0507983,0.9029218,0.8011611];
handles.tdse_run=uicontrol(f,'unit','normalized','BackgroundColor',[-1,-1,-1],'Enable','on','FontAngle','normal','FontName','Tahoma','FontSize',[12],'FontUnits','points','FontWeight','normal','ForegroundColor',[-1,-1,-1],'HorizontalAlignment','center','ListboxTop',[],'Max',[1],'Min',[0],'Position',[0.0537229,0.0406386,0.1046183,0.053701],'Relief','default','SliderStep',[0.01,0.1],'String','Calculate','Style','pushbutton','Value',[0],'VerticalAlignment','middle','Visible','on','Tag','tdse_run','Callback','tdse_run_callback(handles)')


f.visible = "on";


//////////
// Callbacks are defined as below. Please do not delete the comments as it will be used in coming version
//////////

// defining a potential function.
function [y]=g(x)
    y = (x - 4.5) .^ 2
endfunction

// initial conditions.
function [y]=f(x)
    y = exp(- (x - 5) .^ 2)
endfunction

function tdse_run_callback(handles)
//Write your callback for  tdse_run  here
    
    // renews the axes when plotting again.
    delete(handles.tdse_axes.children)
    
    // constants
    hcut = 1
    m = 1
    
    // size of the grid matrices.
    lent = 3000
    lenx = 51
    
    // see what works best.
    delt = 0.001
    
    // The Matrices.
    R = zeros(lent, lenx)
    I = zeros(lent, lenx)
    
    // tailor-made for this application.
    x = linspace(0, 10, lenx)
    v = linspace(0, 10, lenx)
    
    delx = x(2)
    V = g(v)
    
    R(1, :) = f(x)
    I(1, :) = 0
    
    // boundary condition.
    R(:, 1) = 0; R(:, lenx) = 0;
    I(:, 1) = 0; I(:, lenx) = 0;
    
    // solving the PDE.
    for j = 2 : lent
        for i = 2 : lenx - 1
            I(j, i) = I(j - 1, i) - (V(i) .* R(j - 1, i) .* delt ./ hcut) + (hcut .* delt ./ (2 .* m .* delx .* delx)) .* (R(j - 1, i - 1) - 2 .* R(j - 1, i) + R(j - 1, i + 1))
            R(j, i) = R(j - 1, i) + (V(i) .* I(j - 1, i) .* delt ./ hcut) - (hcut .* delt ./ (2 .* m .* delx .* delx)) .* (I(j - 1, i - 1) - 2 .* I(j - 1, i) + I(j - 1, i + 1))
        end
    end
    
    // plotting all iterations.
    //for k = 1 : 20 : lent
    //    plot(x, R(k, :), '-b')
    //    plot(x, I(k, :), ':r')      
    //end
    
    for k = 1 : 10 : lent
        plot(x, (R(k, :) .* R(k, :) + I(k, :) .* I(k, :)), '-b')     
    end
    
    // see the relevant paper for theoretical formulation.

endfunction

