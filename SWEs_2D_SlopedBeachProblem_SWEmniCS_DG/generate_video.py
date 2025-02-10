# Install necessary libraries for video generation
# !pip install imageio
# !pip install imageio[ffmpeg]

import imageio.v2 as imageio  # Import imageio for handling image sequences

def generate_video(output_filename="simulation_wetdry_ramp_weakform_as_swemnics_dg.mp4", fps=12):
    """
    Generates a video from a sequence of simulation frames.

    Parameters:
    - output_filename (str): Name of the output video file.
    - fps (int): Frames per second for the video.

    This function collects all PNG frames stored in the 'frames/' directory, 
    sorts them in order, and compiles them into an MP4 video file.
    """
    print('Running generate_video to see the simulation')

    # Collect all frame files, assuming filenames follow "frame_0000.png", "frame_0001.png", etc.
    frame_files = sorted([f"frames/frame_{i:04d}.png" for i in range(1009)])

    # Create and write the video file
    with imageio.get_writer(output_filename, fps=fps) as writer:
        for frame_file in frame_files:
            image = imageio.imread(frame_file)  # Read image frame
            writer.append_data(image)  # Append frame to video

    print(f"Successfully saved video as {output_filename}")
