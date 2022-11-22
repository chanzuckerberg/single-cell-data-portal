import { IconNames } from "@blueprintjs/icons";
import { SetSearchValueFn } from "src/components/common/Filter/components/FilterSearch/common/useFilterSearch";
import { ViewSearch } from "./style";

interface Props {
  className?: string;
  searchValue: string;
  setSearchValue: SetSearchValueFn;
}

export default function FilterSearch({
  className,
  searchValue,
  setSearchValue,
}: Props): JSX.Element {
  return (
    <ViewSearch
      autoFocus
      className={className}
      leftIcon={IconNames.SEARCH}
      onChange={(changeEvent) => setSearchValue(changeEvent.target.value)}
      placeholder="Search"
      value={searchValue}
    />
  );
}
