import os
import shutil

def clean_pycache(project_root='.'):
    """
    Recursively finds and removes all __pycache__ directories
    starting from the project_root.

    Args:
        project_root (str): The path to the project's root directory.
                            Defaults to the current directory.
    """
    removed_count = 0
    for root, dirs, files in os.walk(project_root, topdown=True):
        if "__pycache__" in dirs:
            pycache_path = os.path.join(root, "__pycache__")
            print(f"Removing: {pycache_path}")
            try:
                shutil.rmtree(pycache_path)
                removed_count += 1
                # Remove __pycache__ from dirs to prevent descending into it after removal
                dirs.remove("__pycache__") 
            except OSError as e:
                print(f"Error removing {pycache_path}: {e}")

    if removed_count > 0:
        print(f"\nSuccessfully removed {removed_count} __pycache__ director(y/ies).")
    else:
        print("\nNo __pycache__ directories found to remove.")

if __name__ == "__main__":
    current_directory = os.getcwd()
    print(f"Starting __pycache__ cleanup in: {current_directory}")
    clean_pycache(current_directory)
    print("Cleanup process finished.")