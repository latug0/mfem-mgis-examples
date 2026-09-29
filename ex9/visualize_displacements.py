#!/usr/bin/env pvpython

"""
Automatic ParaView visualization of displacement fields from a .pvd file.

Designed for:
    ParaView 5.10.x

Features:
    - Opens a PVD time series
    - Detects 2D/3D geometry
    - Handles 2D displacement vectors by creating a 3D vector
    - Warps mesh by displacement
    - Colors by displacement magnitude
    - Shows surface with edges
    - Leaves ParaView in interactive animation mode

Usage:
    pvpython visualize_displacement.py
    pvpython visualize_displacement.py file.pvd
    pvpython visualize_displacement.py file.pvd --warp 5
"""

import sys
import argparse
import math

from paraview.simple import *


# ============================================================
# User configuration
# ============================================================

DEFAULT_FILE = "displacements/displacements.pvd"

DISPLACEMENT_FIELD = "u"

# Default deformation amplification
WARP_SCALE = 1.0


# ============================================================
# Command line
# ============================================================

parser = argparse.ArgumentParser(
    description="Visualize displacement PVD files in ParaView"
)

parser.add_argument(
    "filename",
    nargs="?",
    default=DEFAULT_FILE,
    help="Input .pvd file"
)

parser.add_argument(
    "--warp",
    type=float,
    default=WARP_SCALE,
    help="Warp scale factor"
)

parser.add_argument(
    "--video",
    action="store_true",
    help="Export animation as an MP4"
)

parser.add_argument(
    "--fps",
    type=int,
    default=20,
    help="Frames per second for exported video"
)

args = parser.parse_args()

filename = args.filename
warp_scale = args.warp


# ============================================================
# Load dataset
# ============================================================

print("Opening:", filename)

reader = PVDReader(
    FileName=filename
)
UpdatePipeline(proxy=reader)

## Could lead to errros if field named differently
print("Assuming displacement field:", DISPLACEMENT_FIELD)

# ============================================================
# Detect geometry dimension
# ============================================================
info = reader.GetDataInformation()
bounds = info.GetBounds()

dx = bounds[1] - bounds[0]
dy = bounds[3] - bounds[2]
dz = bounds[5] - bounds[4]

scale = max(dx, dy, dz)

# Geometry is considered 2D if thickness is negligible
is_2D = abs(dz) < 1e-8 * scale


if is_2D:
    print("Detected: 2D geometry")
else:
    print("Detected: 3D geometry")


# ============================================================
# Prepare displacement vector
# ============================================================

if is_2D:

    print("Creating 3D displacement vector")

    calculator = Calculator(
        Input=reader
    )

    calculator.ResultArrayName = "u_3D"

    # Works with standard ParaView component notation
    calculator.Function = (
        "u_X*iHat + u_Y*jHat + 0*kHat"
    )

    UpdatePipeline()

    displacement_array = "u_3D"

    source = calculator

else:

    source = reader

    displacement_array = DISPLACEMENT_FIELD


# ============================================================
# Warp mesh
# ============================================================

warp = WarpByVector(
    Input=source
)

warp.Vectors = [
    "POINTS",
    displacement_array
]

warp.ScaleFactor = warp_scale

UpdatePipeline()


# ============================================================
# Compute displacement magnitude
# ============================================================

calculator_mag = Calculator(
    Input=warp
)

calculator_mag.ResultArrayName = "DisplacementMagnitude"

calculator_mag.Function = (
    "mag(" + displacement_array + ")"
)

UpdatePipeline()


# ============================================================
# Display
# ============================================================

renderView = GetActiveViewOrCreate("RenderView")

display = Show(
    calculator_mag,
    renderView
)


# Surface with edges
display.Representation = "Surface With Edges"


# Color by displacement magnitude
ColorBy(
    display,
    ("POINTS", "DisplacementMagnitude")
)

display.RescaleTransferFunctionToDataRange(
    True,
    False
)


## Show color bar
#display.SetScalarBarVisibility(
#    renderView,
#    True
#)
#
#lut = GetColorTransferFunction("DisplacementMagnitude")
#scalarBar = GetScalarBar(lut, renderView)
#
#scalarBar.Title = "|u|"
#scalarBar.ComponentTitle = ""
#
#scalarBar.Orientation = "Vertical"
#
#scalarBar.WindowLocation = "Any Location"
#
#scalarBar.Position = [0.88, 0.15]
#
#scalarBar.ScalarBarLength = 0.5
#
#scalarBar.TitleFontSize = 16
#scalarBar.LabelFontSize = 14
#
#scalarBar.AutomaticLabelFormat = 0
#scalarBar.LabelFormat = "%.2e"

# ============================================================
# Camera
# ============================================================

if is_2D:

    print("Using 2D camera")

    renderView.CameraPosition = [
        (bounds[0] + bounds[1]) / 2,
        (bounds[2] + bounds[3]) / 2,
        10 * scale
    ]

    renderView.CameraFocalPoint = [
        (bounds[0] + bounds[1]) / 2,
        (bounds[2] + bounds[3]) / 2,
        0
    ]

    renderView.CameraViewUp = [
        0,
        1,
        0
    ]

else:

    print("Using 3D camera")

    renderView.CameraPosition = [
        0.45*bounds[0]+0.55*bounds[1],
        (bounds[2] + bounds[3])/2, #(bounds[2]+ bounds[3])/2 * 1.5,
        bounds[5] * 1.3
    ]

    renderView.CameraFocalPoint = [
        (bounds[0] + bounds[1]) / 2,
        (bounds[2] + bounds[3]) / 2,
        (bounds[4] + bounds[5]) / 2
    ]

    renderView.CameraViewUp = [
        0,
        1,
        0
    ]


ResetCamera()

renderView.Update()


# ============================================================
# Animation setup
# ============================================================

scene = GetAnimationScene()

scene.UpdateAnimationUsingDataTimeSteps()
timesteps = scene.TimeKeeper.TimestepValues
nframes = len(timesteps)

print(f"{nframes} animation frames detected.")

print("Timesteps detected:")
print(scene.TimeKeeper.TimestepValues)


# ============================================================
# Final message
# ============================================================

# ============================================================
# Optional video export
# ============================================================

if args.video:

    video_name = "animation.mp4"
    renderView.ViewSize = [1920, 1080]
    Render()
    print(f"Writing {video_name} ...")
    WriteAnimation(
        "frames/frame.png",
        viewOrLayout=renderView,
        ImageResolution=(1920, 1080),
        )
    print("Done.")
    import shutil
    import subprocess

    if shutil.which("ffmpeg"):
        subprocess.run([
            "ffmpeg",
            "-y",
            "-framerate", str(args.fps),
            "-i", "frames/frame.%04d.png",
            "-c:v", "libx264",
            "-crf", "22",
            "-pix_fmt", "yuv420p",
            video_name,
        ])


print("")
print("Visualization ready.")
print("Press Play in ParaView to animate.")
print("")
