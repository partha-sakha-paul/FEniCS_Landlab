# !pip install imageio
# !pip install imageio[ffmpeg]
import imageio.v2 as imageio

def generate_video(output_filename="simulation_dem_cg.mp4", fps=12):
    print('Creating video to see the simulation')
    # Collect all frame files
    frame_files = sorted([f"frames/frame_{i:04d}.png" for i in range(1009)])

    # Write video file
    with imageio.get_writer(output_filename, fps=fps) as writer:
        for frame_file in frame_files:
            image = imageio.imread(frame_file)
            writer.append_data(image)
    
    print(f"Succesfully Video saved as {output_filename}")
