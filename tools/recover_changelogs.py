import subprocess
import re

# Commits obtained previously (oldest to newest for relevant files)
COMMITS_INFO = [
    ("1acca881844faa763c80a16feb91aa0691690697", "2025-05-09"),
    ("b9bd9711f14c28ae882fd25b5512dae02d599386", "2025-05-09"),
    ("a0a753551acb7113606320f08a98cce4b94f9de6", "2025-05-09"),
    ("623832154111961bf7a5953ee5bbdfbe4b16d788", "2025-05-10"),
    ("9b4c02d958a99d2321ae77accb9e7b23ce756ac0", "2025-05-11"),
    ("c508a1d1774d1764f8e8c9688a62c4d21f46e619", "2025-05-11"),
    ("77d97d9e9f58f50faf83aa8f797a39d5b762f4b4", "2025-05-11"),
    ("a8cf128fdc39fe988f19c269001259056c3e9ffe", "2025-05-11"),
    ("b5c08604dd7f168111c7e193b5a2929fe60b4f49", "2025-05-12"),
    ("9c5777c612794b3ef5e8965887622ba8e35771f7", "2025-05-12"),
    ("fb270a86c92173678675693ef1d5e9bc784bf5b8", "2025-05-13"),
    ("d8dd9eca5d3ee26dfda998d8efa85d9b2490241c", "2025-05-13"),
    ("e87ab7e5387ab82fcd22692a27fb8be136cd0a6e", "2025-05-13"),
    ("4dcabcc4a4794c25e40ed59e04f41f07466c1640", "2025-05-13"),
    ("1524f93e29f0cdacfc85891f0037aad4d1104a3d", "2025-05-13"),
    ("d9a2be501d5d9e387d6aa711750546e5e24987a2", "2025-05-13"),
]

SEARCH_TERM = "2025-05"
# These are pathspecs for git ls-tree, it will find files under these directories.
# The Python filter further ensures only .md files within these exact paths are processed.
PATHSPECS_FOR_GIT_LS_TREE = ["ai/stories", "docs"]

def run_git_command(command):
    try:
        result = subprocess.run(command, capture_output=True, text=True, check=False, shell=False)
        if result.returncode != 0:
            # git ls-tree can return 1 if a pathspec doesn't match, not necessarily a critical error for the script.
            # git show will error if file not in commit.
            # Print stderr for debugging but don't always stop the script.
            if "did not match any file(s) known to git" not in result.stderr:
                 print(f"Warning running command '{' '.join(command)}': {result.stderr}")
            return None # Indicate an issue or no output
        return result.stdout.strip()
    except subprocess.CalledProcessError as e:
        print(f"Error running command '{' '.join(command)}': {e.stderr}")
        return None
    except FileNotFoundError:
        print(f"Error: git command not found. Ensure git is installed and in your PATH.")
        # This is a fatal error for the script, so re-raise or exit
        raise
    except Exception as e:
        print(f"An unexpected error occurred with command '{' '.join(command)}': {e}")
        return None

def main():
    print(f"Searching for '{SEARCH_TERM}' in markdown files within specified paths across recent commits...")

    for commit_hash, commit_date in COMMITS_INFO:
        print(f"\n{'='*60}")
        print(f"Commit: {commit_hash} (Date: {commit_date})")
        print(f"{'='*60}")

        cmd_list_files = ["git", "ls-tree", "-r", "--name-only", commit_hash] + PATHSPECS_FOR_GIT_LS_TREE
        files_output = run_git_command(cmd_list_files)
        
        # If files_output is None (error or no files from ls-tree for that pathspec), 
        # or empty string (no files found under those pathspecs in this commit)
        if not files_output:
            print("  (No .md files found under ai/stories/ or docs/ pathspecs in this commit)")
            continue

        # Filter further to ensure we only process .md files directly under ai/stories or docs
        # and not in sub-sub-directories unless intended by PATHSPECS_FOR_GIT_LS_TREE
        all_files_in_commit_tree_paths = files_output.splitlines()
        relevant_files = [
            f for f in all_files_in_commit_tree_paths
            if f.endswith(".md") and (f.startswith("ai/stories/") or f.startswith("docs/"))
        ]

        if not relevant_files:
            print("  (No relevant .md files matched after filtering in this commit)")
            continue

        commit_had_any_match = False
        for file_path in relevant_files:
            # Ensure file_path is not empty string if splitlines somehow yields one
            if not file_path:
                continue
                
            cmd_show_file = ["git", "show", f"{commit_hash}:{file_path}"]
            file_content = run_git_command(cmd_show_file)

            if file_content is None: # Error showing file (e.g., it was deleted later but ls-tree saw it)
                print(f"  Could not retrieve content for: {file_path}")
                continue

            matches_in_file = []
            for i, line in enumerate(file_content.splitlines()):
                if SEARCH_TERM in line:
                    matches_in_file.append(f"    L{i+1}: {line.strip()}")
            
            if matches_in_file:
                commit_had_any_match = True 
                print(f"\n  >>> File: {file_path}")
                for match_line in matches_in_file:
                    print(match_line)
        
        if not commit_had_any_match:
            print("  (No matches for anTerm in relevant files in this commit)") # Typo: anTerm -> SEARCH_TERM

if __name__ == "__main__":
    main() 