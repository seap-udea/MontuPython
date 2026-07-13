import os
import montu

def format_size(path):
    try:
        size_bytes = os.path.getsize(path)
        if size_bytes >= 1024 * 1024:
            return f"{size_bytes / (1024 * 1024):.2f} MB"
        elif size_bytes >= 1024:
            return f"{size_bytes / 1024:.2f} kB"
        else:
            return f"{size_bytes} bytes"
    except Exception:
        return "unknown"

# Get name of the current stellar catalogue
catalogue_name = montu.stars.STELLAR_CATALOGUE

# Build filenames
visible_filename = catalogue_name.replace('.csv', '_visible.csv')
bright_filename = catalogue_name.replace('.csv', '_bright.csv')

# Get absolute paths to files inside the package data directory
original_path = montu.Util._data_path(catalogue_name)
visible_path = montu.Util._data_path(visible_filename)
bright_path = montu.Util._data_path(bright_filename)

print(f"Loading full stellar catalogue...")
allstars = montu.Stars()
original_count = allstars.number
original_size_str = format_size(original_path)

print(f"Generating visible subset (Vmag in [-3, 6.5])...")
visible = allstars.get_stars(Vmag=[-3, 6.5])
visible.data.to_csv(visible_path, index=False)
visible_count = visible.number
visible_size_str = format_size(visible_path)

print(f"Generating bright subset (Vmag in [-3, 3.5])...")
bright = allstars.get_stars(Vmag=[-3, 3.5])
bright.data.to_csv(bright_path, index=False)
bright_count = bright.number
bright_size_str = format_size(bright_path)

print("\n--- Summary of Stellar Catalogue Subsets ---")
print(f"Original catalogue ({catalogue_name}):")
print(f"  Stars count: {original_count}")
print(f"  File size:   {original_size_str}")

print(f"\nVisible subset ({visible_filename}):")
print(f"  Stars count: {visible_count}")
print(f"  File size:   {visible_size_str}")

print(f"\nBright subset ({bright_filename}):")
print(f"  Stars count: {bright_count}")
print(f"  File size:   {bright_size_str}")
