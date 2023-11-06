# beachball
Leaflet Plugin for Beachballs / Moment Tensors

```js
<!-- After Leaflet -->
<script src="L.Beachball.js"></script>

<script>

const mech = { rr: m.Mrr, tt: m.Mtt, pp: m.Mpp, rt: m.Mrt, rp: m.Mrp, tp: m.Mtp }
const tcolor =  'black'
const bb = L.beachballMarker(L.latLng(m.Latitude, m.Longitude), {
     mech, tcolor
})

</script>

```
