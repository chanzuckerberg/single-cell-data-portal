GPT_GENEGUIDE_DESCRIPTION_SYSTEM_ROLE = "You are a knowledgeable cell biologist and genomics expert that has professional experience writing and curating accurate and informative descriptions of gene ontology terms."


def GPT_GOTERM_DESCRIPTION_USER_ROLE(name: str):
    return f'I am making a knowledge-base about gene ontology (GO) terms. Each GO term from the Gene Ontology and will have its own page with a detailed description of that term. Please write me a description for "{name}". Please return only the description and no other dialogue. The description should include information about biological role of the term \'s function. The description should be at least three paragraphs long.'


def GPT_GOTERM_SEO_DESCRIPTION_USER_ROLE(description: str):
    return f"Here is the description for a gene ontology term's information page:\n\n{description}\n\nPlease write me a short summary of the above description for Search Engine Optimization (SEO). Please return only the summary and no other dialogue. The description should be one-to-two sentences."
