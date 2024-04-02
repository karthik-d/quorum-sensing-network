extensions [ vid sound ]

globals [ Kcell Kpull Kpush Ksubs inv-gamma circle-radius counter counter2 counter3 ]

breed [ neuralcells neuralcell ]
breed [ walls wall ]
neuralcells-own [ FcellX FcellY FpullX FpullY FpushX FpushY FsubsX FsubsY original-heading connectability neurite-count ]
walls-own [ original-heading connectability ]
patches-own [ cluster identity Fperpendicular Fparallel ]
links-own [ cluster1 cluster2 counted real-heading ]

to setup
  clear-all

  set Kcell 5
  set Kpull 1
  set Kpush 5
  set Ksubs adhesion-strength
  set inv-gamma 0.1 ;; 1/gamma damping or friction constant (from PRL paper)
  set circle-radius min list max-pxcor max-pycor

  set-default-shape turtles "circle"

  create-neuralcells cell-number
  [ set color blue set size 1 setxy random-xcor random-ycor
    if circular-area? [ move-to one-of patches with [ (distancexy 0 0) < (circle-radius - 1) ] fd random-float 0.5 ]
    set connectability 0 ]

  ifelse show-cells?
  [ ask neuralcells [ show-turtle ] ]
  [ ask neuralcells [ hide-turtle ] ]

  ask patches
  [ set Fparallel Ksubs
    set Fperpendicular Fparallel
    if nanopattern?
    [ ifelse ((pxcor <= max-pxcor * nanopattern-size / 100) and (pxcor >= min-pxcor * nanopattern-size / 100)) and ((pycor <= max-pycor * nanopattern-size / 100) and (pycor >= min-pycor * nanopattern-size / 100)) ;; align cells inside nanopattern
      [ set Fperpendicular (2 * Fparallel) ]
      [ set Fperpendicular Fparallel ]
    ]
  ]

  ifelse circular-area?
  [ create-walls 2 * ceiling (2 * pi * circle-radius) [ set color yellow ]
    layout-circle walls circle-radius
    __change-topology false false ;; sets wrapping of world OFF
  ]
  [ __change-topology true true ;; sets wrapping of world ON
  ]

  reset-ticks
end


to go
  ask neuralcells
  [ calculate-Fpush
    calculate-Fcell
    calculate-Fpull
    calculate-Fsubs
    calculate-V ]

  ifelse show-cells?
  [ ask neuralcells [ show-turtle ] ]
  [ ask neuralcells [ hide-turtle ] ]

  ifelse show-neurites?
  [ ask links [ show-link ] ]
  [ ask links [ hide-link ] ]

tick

if (ticks = 1000 or ticks = 2000)
[ foreach [60 62 64 65 67 69 71 72] [ [x] -> sound:start-note "XYLOPHONE" x 65 wait 0.2 sound:stop-note "XYLOPHONE" x ]
  stop ]

extend-neurites
end


to calculate-Fcell
  rt random-float 60 - random-float 60
  let Fcell random-float Kcell
  set FcellX (Fcell * dx)
  set FcellY (Fcell * dy)
end


to calculate-Fpull
  ask link-neighbors
  [ set original-heading heading
    face myself rt 180 ]

  if (link-neighbors = nobody) [ set FpullX 0 set FpullY 0 stop ]

  let Mx sum ([dx] of link-neighbors)
  let My sum ([dy] of link-neighbors)

  set FpullX (Kpull * Mx)
  set FpullY (Kpull * My)

  ask link-neighbors [ set heading original-heading ]
end


to calculate-Fpush
  set FpushX 0 set FpushY 0
;  let overlap-cells (turtle-set other neuralcells walls) in-radius size
  let overlap-cells other neuralcells in-radius size ;; comment above line and use this line for much faster simulations (walls don't push cells)

  if (any? overlap-cells)
  [ ask overlap-cells [ face myself ]

    let MMx sum [(size - distance myself) * dx ] of overlap-cells
    let MMy sum [(size - distance myself) * dy ] of overlap-cells

    set FpushX (Kpush * MMx)
    set FpushY (Kpush * MMy)
  ]
end


to calculate-Fsubs
  let FdriveX (FcellX + FpullX)
  let FdriveY (FcellY + FpullY)
  carefully [ set heading (atan FdriveX FdriveY) ] [ ] ;; overall heading for Fdrive

  set FsubsX (Fparallel * dx)
  set FsubsY (Fperpendicular * dy)
end


to calculate-V
  let FdriveX (FcellX + FpullX + FpushX)
  let FdriveY (FcellY + FpullY + FpushY)
  let FresistX (FsubsX)
  let FresistY (FsubsY)

  let FallX (FdriveX - FresistX)
  let FallY (FdriveY - FresistY)

  if ((abs FdriveX - abs FresistX) < 0) [ set FallX 0 ] ;; threshold forces
  if ((abs FdriveY - abs FresistY) < 0) [ set FallY 0 ]

  let Fall sqrt (FallX * FallX + FallY * FallY)

  carefully [ set heading (atan FallX FallY) fd inv-gamma * Fall ] [ ] ;; overall heading
end


to extend-neurites
  ask neuralcells
  [ if ((random-float 100) <= 1 and neurite-count < 10 )
    [ let candidate-cells (turtle-set other neuralcells walls) in-radius 32
      if (candidate-cells = nobody) [ stop ]

      ask candidate-cells
      [ set original-heading heading
        face myself rt 180

        let Phr 1 / (distance myself ^ 2)

        if nanopattern? ;; if nanopattern? and (abs (distance myself * dy) >= 5) [ set Phr 0 ]
        [ ifelse abs ((distance myself) * dy) >= 5
          [ set Phr 0 ]
          [ set Phr 1 / (distance myself ^ 1) ] ]

        if (Phr > random-float 1) [ set connectability 1 ]

        set heading original-heading ]

      let cell-counter min (list (count candidate-cells with [connectability = 1]) 1)
      let connect-cells n-of cell-counter candidate-cells with [connectability = 1]
      create-links-with connect-cells [ set color green ]
      ask candidate-cells [ set connectability 0 ]
      set neurite-count neurite-count + 1
    ]
  ]
end


;;;;;;;;;;;;;;;;;;;;;;;
;;; video recording ;;;
;;;;;;;;;;;;;;;;;;;;;;;

to save-video ;; export a movie of the view
    let _recording-save-file-name "movie.mp4"
    let _recording-time 1000 ;; time (tick) to stop video
    let _frame-skip (_recording-time / 250) ;; how many frames to skip before saving a frame
    vid:start-recorder ;; 25 frames per second, 10 second video
    repeat _recording-time
    [ go
      if (ticks mod _frame-skip) = 0 [ vid:record-view ]
    ]
    vid:save-recording _recording-save-file-name
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; analyze axon and cluster data ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to analyze-data
  get-axon-data
     wait 1
  find-clusters
     wait 1
  get-fasciculation-data
end




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; characterize fascicles ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to get-fasciculation-data

let dummy word "fascicles_" name-of-files
file-open word dummy ".data"

;  ask links ;; with [ link-length >= 5 ]
;  [
;    file-write link-length file-write link-heading file-print ""
;  ]

ask links
  [
    show-link
    let aaa [identity] of end1
    let bbb [identity] of end2
    set cluster1 min list aaa bbb
    set cluster2 max list aaa bbb
    set counted 0
    set real-heading link-heading
  ]

let fascicles links with [ cluster1 != cluster2 and cluster1 != 0 and cluster2 != 0 ]

ask fascicles with [ counted = 0 ]
[
  let dummy-heading link-heading
  let dummy-color one-of base-colors
  let similar-fascicles fascicles with [ cluster1 = [cluster1] of myself and cluster2 = [cluster2] of myself and counted = 0 ]
  ask similar-fascicles
  [
    set counted 1
    set color dummy-color
    let dummy-diff subtract-headings dummy-heading real-heading
    if abs dummy-diff > 90
    [ set real-heading link-heading - 180 ]
    if real-heading < 0
    [ set real-heading real-heading + 360 ]
  ]

  if count similar-fascicles > 0
  [
    file-write count similar-fascicles
    file-write mean [link-length] of similar-fascicles
    file-write mean-of-headings ([real-heading] of similar-fascicles)
    file-print ""
  ]
]

file-close-all

end


to-report mean-of-headings [headings]
  let x-mean mean map sin headings
  let y-mean mean map cos headings
  if x-mean = 0 and y-mean = 0 [ report random 360 ]
  report atan x-mean y-mean
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; extract axon data into file ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to get-axon-data
  let dummy word "axons_" name-of-files
  file-open word dummy ".data"

  ask links ;; with [ link-length >= 5 ]
  [
    file-write link-length file-write link-heading file-print ""
  ]

  file-close-all

end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; cluster characterization ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to find-clusters
  setup-counting
  loop [
    ;; pick a random patch that isn't in a cluster yet
    let seed one-of patches with [cluster = nobody]
    ;; if we can't find one, then we're done!
    if seed = nobody
    [ show-clusters
      set counter counter - 1

      let dummy word "clusters_" name-of-files
      file-open word dummy ".data"

      repeat counter
      [ let dummy-color one-of base-colors - 2
        let cluster-size (count neuralcells-on patches with [identity = counter2])

        ask one-of patches with [identity = counter2]
        [ if (cluster-size >= 10)
          [ set plabel identity ] ]


        ask patches with [identity = counter2]
        [ ifelse (cluster-size >= 10)
          [ set pcolor dummy-color ]
          [ set pcolor black ] ]

        if (cluster-size >= 10)
        [ file-write cluster-size file-print ""
          set counter3 counter3 + 1 ]

        set counter2 counter2 + 1 ]

      file-close-all
      stop ]

    ;; otherwise, make the patch the "leader" of a new cluster
    ;; by assigning itself to its own cluster, then call
    ;; grow-cluster to find the rest of the cluster
    ask seed
    [ set cluster self
      grow-cluster ]
  ]
  display
end


to setup-counting
  ask neuralcells [ hide-turtle ]
  ask links [ hide-link ]

  set counter 1
  set counter2 1
  set counter3 0

  ask patches
  [ ;; use dark colors so the labels are visible
    set pcolor black ;one-of base-colors - 2
    set plabel ""
    set identity 0
    ;; initially, we're in no cluster
    set cluster nobody ]

  ask patches with [count neuralcells in-radius 1 >= 1]
  [ set pcolor red ]
end


;; once all the clusters have been found, this is called
;; to put numeric labels on them so the user can see
;; that the clusters were identified correctly
to show-clusters
  loop
  [ ;; pick a random patch we haven't given an identity number yet
    let p one-of patches with [identity = 0 and pcolor = red]
    if p = nobody
      [ stop ]
    ;; give all patches in the chosen patch's cluster
    ;; the same label
    ask p
    [ ask patches with [cluster = [cluster] of myself]
      [ set identity counter ] ]
    set counter counter + 1 ]
end


to grow-cluster  ;; patch procedure
  ask neighbors4 with [(cluster = nobody) and
;  ask patches in-radius 1.5 with [(cluster = nobody) and
    (pcolor = [pcolor] of myself)]
  [ set cluster [cluster] of myself
    grow-cluster ]
end
@#$#@#$#@
GRAPHICS-WINDOW
230
10
840
621
-1
-1
2.5
1
15
1
1
1
0
0
0
1
-120
120
-120
120
1
1
1
ticks
15.0

SWITCH
12
166
201
199
nanopattern?
nanopattern?
1
1
-1000

SLIDER
13
207
202
240
nanopattern-size
nanopattern-size
0
100
100.0
25
1
%
HORIZONTAL

BUTTON
74
276
137
309
setup
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
17
324
80
357
go
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
12
11
201
44
cell-number
cell-number
100
10000
5000.0
100
1
NIL
HORIZONTAL

SWITCH
12
444
199
477
show-cells?
show-cells?
0
1
-1000

SWITCH
12
404
198
437
show-neurites?
show-neurites?
0
1
-1000

BUTTON
121
324
197
357
save video
save-video
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SWITCH
12
125
200
158
circular-area?
circular-area?
0
1
-1000

MONITOR
909
194
1027
239
number of clusters
counter3
17
1
11

TEXTBOX
860
12
1111
49
To save data after a run (cluster and axon data), input a file name and press 'analyze-data'
11
0.0
1

TEXTBOX
96
334
111
352
or
11
0.0
1

SLIDER
13
58
202
91
adhesion-strength
adhesion-strength
1
5
1.0
1
1
NIL
HORIZONTAL

BUTTON
913
134
1017
167
NIL
analyze-data
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
872
50
1060
110
name-of-files
5000cells_nanopattern
1
0
String

@#$#@#$#@
## WHAT IS IT?

Code to simulate generation of neuronal networks for paper by Onur Kilic et al (2018)

## HOW TO USE IT

1. Adjust 'cell-number' and 'adhesion-strength'
2. Choose whether the simulation should run in a circular area
3. Choose whether there should be a directional bias (due to a nanopatterned surface)
4. Press 'setup'
5. Press 'go' or 'save video' to run simulation ('go' runs regular simulation, 'save video' outputs simulation into a video file 'movie.mp4')
6. After the run, if you want to save the data, press 'analyze-data'
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
