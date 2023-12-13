export function formatCitation(str: string) {
  const [authors, journalYear] = str.split("et al.");
  if (!journalYear) return str;
  const formattedYear = `(${journalYear.slice(-4)})`;
  const journal = journalYear.slice(0, -5);
  return `${authors} et al. ${formattedYear} ${journal}`;
}
