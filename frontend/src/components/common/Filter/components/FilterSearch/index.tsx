import { IconNames } from "@blueprintjs/icons";
import { useRef } from "react";
import { SetSearchValueFn } from "src/components/common/Filter/components/FilterSearch/common/useFilterSearch";
import { ClearButtonIcon, ViewSearch } from "./style";

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
  const inputRef = useRef<HTMLInputElement>(null);

  /**
   * Clears search value and sets focus on the search input element.
   */
  const onClearSearch = (): void => {
    setSearchValue("");
    inputRef?.current?.focus();
  };

  return (
    <ViewSearch
      autoFocus
      className={className}
      inputRef={inputRef}
      leftIcon={IconNames.SEARCH}
      onChange={(changeEvent) => setSearchValue(changeEvent.target.value)}
      placeholder="Search"
      rightElement={
        searchValue ? (
          <ClearButtonIcon
            onClick={onClearSearch}
            sdsSize="small"
            sdsType="secondary"
            icon="XMarkCircle"
          />
        ) : undefined
      }
      value={searchValue}
    />
  );
}
