GPT_CELLTYPE_DESCRIPTION_SYSTEM_ROLE = "You are a knowledgeable cell biologist that has professional experience writing and curating accurate and informative descriptions of cell types."


def GPT_CELLTYPE_DESCRIPTION_USER_ROLE(cname):
    return f'I am making a knowledge-base about cell types. Each cell type is a term from the Cell Ontology and will have its own page with a detailed description of that cell type and its function. Please write me a description for "{cname}". Please return only the description and no other dialogue. The description should include information about the cell type\'s function. The description should be at least three paragraphs long.'
