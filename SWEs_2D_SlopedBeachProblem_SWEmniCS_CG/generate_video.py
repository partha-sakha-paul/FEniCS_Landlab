# !pip install imageio
# !pip install imageio[ffmpeg]
import imageio.v2 as imageio

def generate_video(output_filename="simulation_wetdry_ramp_weakform_as_swemnics_cg.mp4", fps=12):
    """
    Generate a video from simulation frames.

    Parameters:
    - output_filename (str): Name of the output video file.
    - fps (int): Frames per second for the output video.

    Returns:
    - None (Saves the video file in the working directory).
    """
    print('Running generate_video to see the simulation')

    # Collect all frame image files, assuming they are named as 'frame_0000.png', 'frame_0001.png', ...
    frame_files = sorted([f"frames/frame_{i:04d}.png" for i in range(1009)])

    # Create and write frames to the video file
    with imageio.get_writer(output_filename, fps=fps) as writer:
        for frame_file in frame_files:
            image = imageio.imread(frame_file)  # Read each frame image
            writer.append_data(image)  # Append to the video

    print(f"Successfully saved video as {output_filename}")
