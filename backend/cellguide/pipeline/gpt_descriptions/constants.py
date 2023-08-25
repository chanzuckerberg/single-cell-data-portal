GPT_CELLTYPE_DESCRIPTION_SYSTEM_ROLE = "You are a knowledgeable cell biologist that has professional experience writing and curating accurate and informative descriptions of cell types."


def GPT_CELLTYPE_DESCRIPTION_USER_ROLE(name: str):
    return f'I am making a knowledge-base about cell types. Each cell type is a term from the Cell Ontology and will have its own page with a detailed description of that cell type and its function. Please write me a description for "{name}". Please return only the description and no other dialogue. The description should include information about the cell type\'s function. The description should be at least three paragraphs long.'


def GPT_CELLTYPE_SEO_DESCRIPTION_USER_ROLE(description: str):
    return f"Here is the description for a cell type's information page:\n\n{description}\n\nPlease write me a short summary of the above description for Search Engine Optimization (SEO). Please return only the summary and no other dialogue. The description should be one-to-two sentences."
