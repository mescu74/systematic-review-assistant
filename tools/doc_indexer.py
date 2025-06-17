import os
import re
from collections import defaultdict

# --- CONFIGURATION ---

DOCS_ROOT = 'docs'
INDEX_PATH = os.path.join(DOCS_ROOT, 'index.md')

# Define categories and rules for assigning files.
# The tuple is (Category Title, Sub-category Title)
# More specific paths should come first.
CATEGORIES: dict[str, tuple[str, str]] = {
    # Specific file assignments
    'architecture.md': ("🏛️ System Architecture", "Overviews"),
    'api-reference.md': ("🏛️ System Architecture", "API & Service Specifications"),
    'architecture/data-models.md': ("🏛️ System Architecture", "Architecture Components"),
    'architecture/database-schema.md': ("🏛️ System Architecture", "Architecture Components"),
    'architecture/unified-project-structure.md': ("🏛️ System Architecture", "Architecture Components"),
    'architecture/tech-stack.md': ("🏛️ System Architecture", "Architecture Components"),
    'architecture/coding-standards.md': ("💻 Development & Operations", "Core Development Guides"),
    'architecture/testing-strategy.md': ("💻 Development & Operations", "Core Development Guides"),
    'guides/prompt-engineering-guide.md': ("💻 Development & Operations", "Core Development Guides"),
    'guides/streamlit-apptest-framework-guide.md': ("💻 Development & Operations", "Core Development Guides"),
    'prd.md': ("📋 Project & Product Management", "Product Requirements Documents (PRDs)"),
    'prd-benchmark-may.md': ("📋 Project & Product Management", "Product Requirements Documents (PRDs)"),
    'prd-recovery.md': ("📋 Project & Product Management", "Product Requirements Documents (PRDs)"),
    'prd-resolver.md': ("📋 Project & Product Management", "Product Requirements Documents (PRDs)"),
    'project-brief.md': ("📋 Project & Product Management", "Project Briefs & Overviews"),
    'project-brief.recovery.md': ("📋 Project & Product Management", "Project Briefs & Overviews"),
    'sr-metrics.md': ("💻 Development & Operations", "Benchmarking & Metrics"),
    'benchmark/benchmark-protocol.md': ("💻 Development & Operations", "Benchmarking & Metrics"),
    'benchmark/benchmark-workflows.md': ("💻 Development & Operations", "Benchmarking & Metrics"),
    'benchmark/o3-xml-criteria-homelessness.md': ("💻 Development & Operations", "Benchmarking & Metrics"),
    'environment-vars.md': ("💻 Development & Operations", "Core Development Guides"),
    'naming-conventions.md': ("💻 Development & Operations", "Core Development Guides"),

    # Directory-based assignments
    'architecture': ("🏛️ System Architecture", "Architecture Components"),
    'guides': ("💻 Development & Operations", "Core Development Guides"),
    'prd': ("📋 Project & Product Management", "Product Requirements Documents (PRDs)"),
    'stories': ("📋 Project & Product Management", "User Stories"),
    'benchmark': ("💻 Development & Operations", "Benchmarking & Metrics"),
    'chats': ("📜 Archives & Miscellaneous", "Archived Chat Logs"),
}
DEFAULT_CATEGORY = ("📜 Archives & Miscellaneous", "Other Documents")

# --- HELPER FUNCTIONS ---

def get_doc_title(content: str, filepath: str) -> str:
    """Extracts title from H1, or generates a fallback title."""
    match = re.search(r"^#\s+(.*)", content, re.MULTILINE)
    if match:
        return match.group(1).strip().replace('**', '')
    name = os.path.basename(filepath).replace('.md', '').replace('-', ' ').replace('_', ' ')
    return name.title()

def generate_description(content: str, title: str) -> str:
    """Generates a brief, one-sentence description from the document's content."""
    content_no_headings = re.sub(r"^[#`].*", "", content, flags=re.MULTILINE).strip()
    content_no_headings = re.sub(r"---.*---", "", content_no_headings, flags=re.DOTALL)
    
    first_paragraph = content_no_headings.strip().split('\n\n')[0].replace('\n', ' ')
    if first_paragraph:
        # Clean up markdown and limit length
        clean_text = re.sub(r'[\\*#_`]', '', first_paragraph)
        words = clean_text.split()
        desc = ' '.join(words[:25])
        return desc + ('...' if len(words) > 25 else '.')
    return f"This document provides details on {title.lower()}."

def parse_existing_index(index_path: str) -> dict[str, tuple[str, str]]:
    """Parses existing index.md to preserve titles and descriptions."""
    if not os.path.exists(index_path):
        return {}
    
    with open(index_path, 'r', encoding='utf-8') as f:
        content = f.read()

    indexed_files = {}
    pattern = re.compile(r"\[(?P<title>.*?)]\((?P<path>.*?)\)\*?\*?:?\s*(?P<desc>.*?)(?=\n-|\n\n##|\Z)", re.DOTALL)
    
    for match in pattern.finditer(content):
        path = match.group('path')
        title = match.group('title').replace('**', '')
        desc = match.group('desc').strip().replace('\n', ' ')
        if path and title and desc:
            indexed_files[path] = (title, desc)
    return indexed_files

def get_category(filepath: str) -> tuple[str, str]:
    """Assigns a file to a category based on the configuration."""
    rel_path = os.path.relpath(filepath, DOCS_ROOT)
    
    # Check for full path matches first
    if rel_path in CATEGORIES:
        return CATEGORIES[rel_path]

    # Check for directory/prefix matches
    parts = rel_path.split(os.sep)
    if parts[0] in CATEGORIES:
        return CATEGORIES[parts[0]]
        
    return DEFAULT_CATEGORY

# --- MAIN LOGIC ---

def main():
    print("Starting documentation indexing...")
    existing_index = parse_existing_index(INDEX_PATH)
    
    # 1. Find all markdown files
    all_md_files = []
    for root, _, files in os.walk(DOCS_ROOT):
        for file in files:
            if file.endswith('.md') and not file.endswith('.bak') and os.path.basename(file) != 'index.md':
                all_md_files.append(os.path.join(root, file))

    # 2. Categorize all documents
    categorized_docs: dict[str, dict[str, list[dict[str, str]]]] = defaultdict(lambda: defaultdict(list))
    for filepath in all_md_files:
        cat, sub_cat = get_category(filepath)
        relative_path = './' + os.path.relpath(filepath, DOCS_ROOT)

        if relative_path in existing_index:
            title, desc = existing_index[relative_path]
        else:
            try:
                with open(filepath, 'r', encoding='utf-8') as f:
                    content = f.read()
                title = get_doc_title(content, filepath)
                desc = generate_description(content, title)
            except Exception:
                title = os.path.basename(filepath)
                desc = "Could not read or parse file."

        categorized_docs[cat][sub_cat].append({'title': title, 'path': relative_path, 'desc': desc})
    
    # 3. Build the new index content
    output_content = [
        "# Documentation Index",
        "This document serves as the central hub for all project documentation. It is organized to help you quickly find information related to project management, architecture, development practices, and historical context."
    ]

    sorted_categories = sorted(categorized_docs.items())

    for category, sub_categories in sorted_categories:
        output_content.append(f"\n---\n\n## {category}")
        sorted_sub_categories = sorted(sub_categories.items())
        
        for sub_category, docs in sorted_sub_categories:
            output_content.append(f"\n### {sub_category}")
            
            # Sort documents alphabetically by title for consistency
            sorted_docs = sorted(docs, key=lambda x: x['title'])
            for doc in sorted_docs:
                # Using a list format for a cleaner look
                output_content.append(f"- **[{doc['title']}]({doc['path']})**: {doc['desc']}")

    # 4. Write the new index file
    with open(INDEX_PATH, 'w', encoding='utf-8') as f:
        f.write('\n'.join(output_content))

    print(f"Successfully updated '{INDEX_PATH}'")
    print(f"Indexed {len(all_md_files)} documents into {len(categorized_docs)} categories.")

if __name__ == "__main__":
    main() 