import Link from "next/link";
import { ROUTES } from "src/common/constants/routes";
import Step from "./components/Step";
import { ColumnOne, ColumnTwo, Details, Header, StepThree, StepTwo, StyledStepOne, StyledStepThree, StyledStepTwo, Wrapper} from "./style";

interface Props {
  tissueSelected: boolean;
  isLoading: boolean;
  geneSelected: boolean;
}

let tissueHasLoadedOnce = false;
let geneHasLoadedOnce = false;

export default function GetStarted({ tissueSelected, isLoading, geneSelected }: Props): JSX.Element {
  if(!tissueHasLoadedOnce && tissueSelected && !isLoading){
    tissueHasLoadedOnce = true;
  }
  if(!geneHasLoadedOnce && geneSelected && !isLoading){
    geneHasLoadedOnce = true;
  }
  return (
    <Wrapper style={ tissueHasLoadedOnce && geneHasLoadedOnce ? {display: "none"} : {}}>

      <ColumnOne style={ tissueHasLoadedOnce ? { visibility: "hidden"} : {} }>
        <StyledStepOne>
          <Step step={1} details="Add Tissues" />
        </StyledStepOne>
      </ColumnOne>

      <ColumnTwo style={ geneHasLoadedOnce ? { visibility: "hidden"} : {} }>
        <StyledStepTwo>
          <Step step={2} details="Add Genes" />
        </StyledStepTwo>
        <StyledStepThree>
          <Step step={3} details="Explore Gene Expression" />
        </StyledStepThree>
      </ColumnTwo>


      {/* <Header>Getting Started</Header> */}
      {/* <Details>
        Use the Add Tissue and Add Gene buttons to find where genes are
        expressed, powered by data from the{" "}
        <Link href={ROUTES.COLLECTIONS} passHref>
          <a href="passHref" rel="noopener" target="_blank">
            CELLxGENE Discover
          </a>
        </Link>
        .
      </Details> */}
    </Wrapper>
  );
}
