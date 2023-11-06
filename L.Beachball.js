
import { plotmt } from './beachball.js'

L.BeachballIcon = L.Icon.extend({
    options: {

    },

    collect: [],
    colors: [],
    size: 200,

    createIcon: function() {
        var e = document.createElement("canvas");
        this._setIconStyles(e, "icon");
        var size = this.options.iconSize;
        const out = plotmt(this.options.mech, {
            size: size[0],
            tcolor: this.options.tcolor,
            pcolor: this.options.pcolor
        })
        this.collect = out.collect;
        this.colors = out.colors;
        e.width = size[0];
        e.height = size[1];
        this.ctx = e.getContext("2d");
        this.draw(this.ctx, size[0], size[1]);

        return e;
    },

    draw: function(ctx, w, h) {
        if(!ctx) return;
        ctx.clearRect(0, 0, w, h);
        ctx.save()
        ctx.translate(w/2, h/2);
        ctx.beginPath();
        for(let j = 0; j < this.collect.length; j++) {
            let x = this.collect[j].x
            let y = this.collect[j].y
            ctx.beginPath()
            ctx.fillStyle = this.colors[j]
            ctx.strokeStyle = this.colors[j]//'rgba(0,0,0,1.0)'
            ctx.lineWidth = 0
            ctx.moveTo(x[0], -y[0] )
            for(let i = 1; i < x.length; i++) {
                ctx.lineTo(x[i], -y[i])
            }
            ctx.closePath()
            ctx.fill()
            ctx.stroke()
        }
    },
});

L.BeachballMarker = L.Marker.extend({
    options: {},
    initialize: function(pos, options) {
        this.setLatLng(pos)
        L.setOptions(this, options)
    },
    pos: function () {
        return this.getLatLng();
    },
});

L.beachballMarker = function(pos, options) {
    const iconSizeDefault = 20;
    options.tcolor = ("tcolor" in options) ? options.tcolor : 'black'
    options.pcolor = ("pcolor" in options) ? options.pcolor : 'white'
    options.iconSize = ("iconSize" in options) ? options.iconSize : [iconSizeDefault,iconSizeDefault]
    options.mech = ("mech" in options) ? options.mech : {xx:0,yy:0,zz:0,xy:0,xz:0,yz:1}

    options.iconSize = ("size" in options) ? options.size : [iconSizeDefault,iconSizeDefault]
    options.iconAnchor = ("iconAnchor" in options) ? options.iconAnchor : [options.iconSize[0]/2, options.iconSize[1]/2]

    options.icon = new L.BeachballIcon(options);

    return new L.BeachballMarker(pos, options)
};
