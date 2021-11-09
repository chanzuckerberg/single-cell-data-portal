import { MenuItem } from "@blueprintjs/core";
import { IItemRendererProps, MultiSelect } from "@blueprintjs/select";
import { useState } from "react";
import { GENES } from "../../mocks/brain";

export default function GeneSearchBar(): JSX.Element {
  const [selectedGenes, setSelectedGenes] = useState([]);

  return (
    <MultiSelect
      onItemSelect={handleItemSelect}
      onRemove={handleItemRemove}
      items={GENES}
      itemRenderer={renderItem}
      tagRenderer={TagRenderer}
      itemsEqual={areGenesEqual}
      selectedItems={selectedGenes}
    />
  );

  function handleItemSelect(gene) {
    if (selectedGenes.includes(gene)) {
      handleItemRemove(gene);
    } else {
      setSelectedGenes([...selectedGenes, gene]);
    }
  }

  function handleItemRemove(gene) {
    setSelectedGenes(
      selectedGenes.filter((selectedGene) => selectedGene !== gene)
    );
  }

  function renderItem(gene, itemRendererProps: IItemRendererProps) {
    return (
      <ItemRenderer
        item={gene}
        isSelected={isGeneSelected(gene)}
        {...itemRendererProps}
      />
    );
  }

  function isGeneSelected(gene) {
    return selectedGenes.includes(gene);
  }
}

function ItemRenderer({
  item,
  handleClick,
  modifiers,
  query,
  isSelected,
}): JSX.Element | null {
  if (!modifiers.matchesPredicate) {
    return null;
  }

  const { name } = item;

  return (
    <MenuItem
      active={isSelected}
      key={name}
      onClick={handleClick}
      text={highlightText(name, query)}
    />
  );
}

function highlightText(text: string, query: string) {
  let lastIndex = 0;
  const words = query
    .split(/\s+/)
    .filter((word) => word.length > 0)
    .map(escapeRegExpChars);
  if (words.length === 0) {
    return [text];
  }
  const regexp = new RegExp(words.join("|"), "gi");
  const tokens: React.ReactNode[] = [];
  while (true) {
    const match = regexp.exec(text);
    if (!match) {
      break;
    }
    const length = match[0].length;
    const before = text.slice(lastIndex, regexp.lastIndex - length);
    if (before.length > 0) {
      tokens.push(before);
    }
    lastIndex = regexp.lastIndex;
    tokens.push(<strong key={lastIndex}>{match[0]}</strong>);
  }
  const rest = text.slice(lastIndex);
  if (rest.length > 0) {
    tokens.push(rest);
  }
  return tokens;
}

function escapeRegExpChars(text: string) {
  return text.replace(/([.*+?^=!:${}()|[\]/\\])/g, "\\$1");
}

function TagRenderer({ name }) {
  return name;
}

function areGenesEqual(geneA, geneB) {
  // Compare only the names (ignoring case) just for simplicity.
  return geneA.name.toLowerCase() === geneB.name.toLowerCase();
}
