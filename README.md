# beachball
Leaflet Plugin for Beachballs / Moment Tensors

```
<!-- After Leaflet -->
<script src="beachball.js"></script>

<script>

const mech = { xx: m.Mrr, yy: m.Mtt, zz: m.Mpp, xy: m.Mrt, xz: m.Mrp, yz: m.Mtp }
const tcolor =  'black'
const bb = L.beachballMarker(L.latLng(m.Latitude, m.Longitude), {
     mech, tcolor
})

</script>

```
