import { useEffect, useState } from "react";
import Step from "./components/Step";
import {
  ColumnOne,
  ColumnTwo,
  StyledStepOne,
  StyledStepThree,
  StyledStepTwo,
  Wrapper,
} from "./style";

interface Props {
  tissueSelected: boolean;
  isLoading: boolean;
  geneSelected: boolean;
}

let tissueHasLoadedOnce = false;
let geneHasLoadedOnce = false;

export default function GetStarted({
  tissueSelected,
  isLoading,
  geneSelected,
}: Props): JSX.Element | null {
  if (!tissueHasLoadedOnce && tissueSelected && !isLoading) {
    tissueHasLoadedOnce = true;
  }
  if (!geneHasLoadedOnce && geneSelected) {
    geneHasLoadedOnce = true;
  }

  const [isClient, setIsClient] = useState(false);

  useEffect(() => {
    setIsClient(true);
  }, []);
  if (isClient)
    return (
      isClient && (
        <Wrapper isHidden={tissueHasLoadedOnce && geneHasLoadedOnce}>
          <ColumnOne>
            <StyledStepOne isHidden={tissueHasLoadedOnce}>
              <Step step={1} details="Add Tissues" />
            </StyledStepOne>
          </ColumnOne>

          <ColumnTwo>
            <StyledStepTwo isHidden={geneHasLoadedOnce}>
              <Step step={2} details="Add Genes" />
            </StyledStepTwo>
            <StyledStepThree
              isHidden={geneHasLoadedOnce && tissueHasLoadedOnce}
            >
              <Step step={3} details="Explore Gene Expression" />
            </StyledStepThree>
          </ColumnTwo>
        </Wrapper>
      )
    );
  return null;
}
