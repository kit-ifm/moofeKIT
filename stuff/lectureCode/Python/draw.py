import matplotlib.pyplot as plt


def draw(ed):

    

#x=ed(:,1:2:size(ed,2))';
#y=ed(:,2:2:size(ed,2))';  

#plot(x,y,pprops{:})

    x = ed[:, 0:ed.shape[1]:2].T
    y = ed[:, 1:ed.shape[1]:2].T
    feplot = plt.plot(x, y)
    plt.setp(feplot, 'color', 'r', 'linewidth', 1)
    plt.axis('equal')
    plt.show
